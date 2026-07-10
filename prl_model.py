import numpy as np

class PLRModel:
    """
    Nonlinear Longtin-Milton DDE model of the Pupillary Light Reflex.

    Governing equation:
        tau_iris * dA/dt = -A(t) + A_star + gamma * [ f(Phi(t - delta)) - f(A_star) ]

    where Phi(t - delta) = A(t - delta) * (1 + s(t - delta)) is the delayed retinal
    light flux, and f(Phi) is the Hill function representing the negative feedback:
        f(Phi) = c * theta^n / (theta^n + Phi^n)

    The scaling factor gamma is chosen such that the linearized local loop gain
    at equilibrium is exactly G:
        gamma = G / |f'(A_star)|

    This ensures the model mathematically conforms to the Lambert W linear theory
    for small perturbations, while saturating physically during large fluctuations.
    """

    def __init__(self, G, proc_noise_sigma=0.0, obs_noise_sigma=0.0, dt=0.001):
        self.dt          = dt
        self.delta       = 0.300
        self.tau_iris    = 0.311
        
        # Hardware Mesopic Lock constraint (aligned to Section 2.1)
        self.A_star      = 12.0

        self.G                = G
        self.proc_noise_sigma = proc_noise_sigma  # Internal biological noise
        self.obs_noise_sigma  = obs_noise_sigma   # CCD camera measurement noise

        # Longtin-Milton Canonical Parameters (aligned to Section 2.6)
        self.n     = 2.0
        self.theta = 12.0
        self.c     = 420.0

        # Derivative of f at equilibrium
        # f'(A) = - c * n * theta^n * A^(n-1) / (theta^n + A^n)^2
        f_prime = -self.c * self.n * (self.theta**self.n) * (self.A_star**(self.n - 1)) / \
                  ((self.theta**self.n + self.A_star**self.n)**2)

        self.gamma = self.G / abs(f_prime)
        self.f_eq  = self._hill_function(self.A_star)

        self.N_delay         = int(round(self.delta / self.dt))
        self.history_buffer  = np.full(self.N_delay, self.A_star)
        self.stimulus_buffer = np.zeros(self.N_delay)
        self.buffer_idx      = 0

        self.A_current = self.A_star

    def _hill_function(self, Phi):
        return self.c * (self.theta**self.n) / (self.theta**self.n + Phi**self.n)

    def reset(self):
        self.history_buffer.fill(self.A_star)
        self.stimulus_buffer.fill(0.0)
        self.buffer_idx = 0
        self.A_current  = self.A_star

    def step(self, stimulus_val):
        A_delayed = self.history_buffer[self.buffer_idx]
        s_delayed = self.stimulus_buffer[self.buffer_idx]

        # 1. PROCESS NOISE (intrinsic mechanical/neural variance before feedback)
        if self.proc_noise_sigma > 0:
            A_delayed += self.proc_noise_sigma * self.A_star * np.random.randn()
            A_delayed = max(A_delayed, 0.1)

        Phi_delayed = A_delayed * (1.0 + s_delayed)

        feedback = self.gamma * (self._hill_function(Phi_delayed) - self.f_eq)

        dA_dt = (1.0 / self.tau_iris) * (-self.A_current + self.A_star + feedback)

        self.A_current += dA_dt * self.dt
        self.A_current = max(self.A_current, 0.1)

        self.history_buffer[self.buffer_idx]  = self.A_current
        self.stimulus_buffer[self.buffer_idx] = stimulus_val
        self.buffer_idx = (self.buffer_idx + 1) % self.N_delay

        # Note: observational noise is not added to the internal state A_current.
        return self.A_current

    def simulate(self, stimulus_array):
        areas = np.zeros(len(stimulus_array))
        for i, stim in enumerate(stimulus_array):
            areas[i] = self.step(stim)
            
        # 2. OBSERVATIONAL NOISE (CCD camera noise on the final output trace)
        if self.obs_noise_sigma > 0:
            noise_array = self.obs_noise_sigma * self.A_star * np.random.randn(len(areas))
            areas += noise_array
            areas = np.maximum(areas, 0.1)  # Physical clamp to prevent negative area readings
            
        return areas