import numpy as np

class PLRModel:
    """
    Nonlinear Delay-Differential Equation (DDE) model of the Pupillary Light Reflex.
    Implements the explicit Euler integrator with a history buffer.
    """
    
    def __init__(self, G, noise_sigma=0.0, dt=0.001):
        # Simulation parameters 
        self.dt = dt
        self.delta = 0.300           # Neural conduction delay (s)
        self.tau_iris = 0.311        # Iris time constant (s)
        self.n = 3                   # Hill coefficient
        self.theta = 50.0            # Half-saturation constant (mm^2)
        self.c = 200.0               # Maximum sphincter area (mm^2)
        self.A_star = 15.0           # Resting pupil area (mm^2)
        
        self.G = G
        self.noise_sigma = noise_sigma
        
        # History buffer setup
        self.N_delay = int(self.delta / self.dt)
        self.history_buffer = np.full(self.N_delay, self.A_star)
        self.buffer_idx = 0
        
        # Compute equilibrium feedback and gain scaling
        self.f_eq = self._hill_function(self.A_star)
        
        # Derivative of Hill function at A_star
        # f'(A) = -c * n * theta^n * A^(n-1) / (theta^n + A^n)^2
        numerator = -self.c * self.n * (self.theta ** self.n) * (self.A_star ** (self.n - 1))
        denominator = (self.theta ** self.n + self.A_star ** self.n) ** 2
        self.f_prime_A_star = numerator / denominator
        
        # Scale optical coupling coefficient alpha to achieve target G
        # G = alpha * |f'(A*)|  =>  alpha = G / |f'(A*)|
        self.alpha = self.G / abs(self.f_prime_A_star)
        
        # State variable
        self.A_current = self.A_star

    def _hill_function(self, A):
        """Nonlinear sphincter response function."""
        # Prevent negative values due to noise from causing numerical issues
        A = max(0.0, A)
        return (self.c * (self.theta ** self.n)) / ((self.theta ** self.n) + (A ** self.n))

    def reset(self):
        """Resets the model to equilibrium state."""
        self.history_buffer.fill(self.A_star)
        self.buffer_idx = 0
        self.A_current = self.A_star

    def step(self, stimulus_val):
        """
        Advances the DDE by one time step dt.
        stimulus_val: proportional increase in effective retinal flux.
        """
        # Retrieve delayed area from history buffer
        A_delayed = self.history_buffer[self.buffer_idx]
        
        # Apply Gaussian additive noise to the delayed area
        if self.noise_sigma > 0:
            noise = self.noise_sigma * self.A_star * np.random.randn()
            A_delayed_noisy = A_delayed + noise
        else:
            A_delayed_noisy = A_delayed
            
        # The light stimulus acts as a multiplier on the effective area 'seen' by the retina
        A_eff = A_delayed_noisy * (1.0 + stimulus_val)
        
        # Compute the nonlinear feedback deviation scaled by alpha
        feedback = self.alpha * (self._hill_function(A_eff) - self.f_eq)
        
        # Explicit Euler integration step
        dA_dt = (1.0 / self.tau_iris) * (-self.A_current + self.A_star + feedback)
        self.A_current += dA_dt * self.dt
        
        # Update history buffer
        self.history_buffer[self.buffer_idx] = self.A_current
        self.buffer_idx = (self.buffer_idx + 1) % self.N_delay
        
        return self.A_current

    def simulate(self, stimulus_array):
        """
        Simulates the model over a provided array of stimulus values.
        Returns an array of pupil areas of the same length.
        """
        areas = np.zeros(len(stimulus_array))
        for i, stim in enumerate(stimulus_array):
            areas[i] = self.step(stim)
        return areas