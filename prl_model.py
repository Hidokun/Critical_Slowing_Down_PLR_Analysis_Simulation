import numpy as np

class PLRModel:
    """
    Pupillary Light Reflex model based on the linearised Longtin-Milton DDE.

    Governing equation (deviation form):
        tau_iris * da/dt = -a(t) + G * a(t - delta) + u(t)

    where a(t) = A(t) - A_star is the deviation from resting pupil area,
    G is the dimensionless loop gain, delta is the neural conduction delay,
    and u(t) is the stimulus input (negative = constriction).

    Characteristic equation: tau_iris * lambda + 1 - G * exp(-lambda * delta) = 0
    Bifurcation at G = 1 (lambda = 0 is a solution).
    Critical slowing down: tau_return = tau_iris / (1 - G) as G -> 1.

    The Hopf bifurcation frequency at G = 1:
        cos(omega_c * delta) = -1  =>  omega_c = pi / delta
        f_c = 1 / (2 * delta) ~= 1.67 Hz for delta = 0.300 s

    Stimulus implementation:
        A light pulse of amplitude S displaces a(t) by -S * A_star
        (constriction: area decreases below A_star during pulse).
    """

    def __init__(self, G, noise_sigma=0.0, dt=0.001):
        # Fixed physiological parameters (Longtin & Milton 1989)
        self.dt         = dt
        self.delta      = 0.300    # Neural conduction delay (s)
        self.tau_iris   = 0.311    # Iris time constant (s)
        self.A_star     = 15.0     # Resting pupil area (mm^2)

        # Gain and noise
        self.G           = G
        self.noise_sigma = noise_sigma

        # History buffer for delayed state
        self.N_delay        = int(round(self.delta / self.dt))
        self.history_buffer = np.zeros(self.N_delay)   # deviations, init at 0
        self.buffer_idx     = 0

        # State variable: deviation from equilibrium
        self.a_current = 0.0

    @property
    def A_current(self):
        """Absolute pupil area (mm^2)."""
        return self.A_star + self.a_current

    def reset(self):
        """Resets model to equilibrium."""
        self.history_buffer.fill(0.0)
        self.buffer_idx = 0
        self.a_current  = 0.0

    def step(self, stimulus_val):
        """
        Advances DDE by one time step dt.

        stimulus_val: float
            Fractional constriction drive. During a pulse, stimulus_val > 0
            produces a negative displacement: u = -stimulus_val * A_star.
            Between pulses, stimulus_val = 0 and u = 0.

        Returns absolute pupil area A(t) = A_star + a(t).
        """
        # Retrieve delayed deviation from history buffer
        a_delayed = self.history_buffer[self.buffer_idx]

        # Additive Gaussian noise on delayed state
        if self.noise_sigma > 0:
            a_delayed += self.noise_sigma * self.A_star * np.random.randn()

        # Stimulus input: light pulse constricts pupil (negative displacement)
        u = -stimulus_val * self.A_star

        # Linearised DDE: tau * da/dt = -a(t) + G * a(t-delta) + u(t)
        da_dt = (1.0 / self.tau_iris) * (
            -self.a_current
            + self.G * a_delayed
            + u
        )

        self.a_current += da_dt * self.dt

        # Write current deviation into history buffer
        self.history_buffer[self.buffer_idx] = self.a_current
        self.buffer_idx = (self.buffer_idx + 1) % self.N_delay

        return self.A_current

    def simulate(self, stimulus_array):
        """
        Simulates model over a stimulus array.
        Returns absolute pupil area array of the same length.
        """
        areas = np.zeros(len(stimulus_array))
        for i, stim in enumerate(stimulus_array):
            areas[i] = self.step(stim)
        return areas