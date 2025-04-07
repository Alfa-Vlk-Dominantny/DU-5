# DU-5
Na ODR- ky


#Uloha 1

import numpy as np
import matplotlib.pyplot as plt

def EulerMethod(f, y0, t, args=()):
    """
    f je tu funkcia pravej strany a y0 počiatocna podmienka a args dalsie podmienky f ktore vlozime do tejto funkcie
    """
    
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(n-1):
        h = t[i+1] - t[i]
        y[i+1] = y[i] + h * f(y[i], t[i], *args)
    return y

def RungeKutta4Order(f, y0, t, args=()):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(n-1):
        h = t[i+1] - t[i]
        k1 = f(y[i], t[i], *args)
        k2 = f(y[i] + k1*h/2, t[i] + h/2, *args)
        k3 = f(y[i] + k2*h/2, t[i] + h/2, *args)
        k4 = f(y[i] + k3*h, t[i] + h, *args)
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4)*h/6
    return y

def pendulum(y, t, eta, omega, A, omega_b):
    theta, omega_t = y
    dydt = [
        omega_t,
        -eta*omega_t - omega**2*np.sin(theta) + A*np.sin(omega_b*t)
    ]

    """
      Zavedena naša rovnica, y tam bude stavovy vektor, kde omega_t je d theta / dt. Tato funkcia vrati teda derivacie y v nejakom case
    """
    return np.array(dydt)


eta = 0.5
omega = 1
A = 1
omega_b_values = [0.2, 0.4, 0.6]
t_max = 500
dt_values = [0.01, 0.005]
y0 = [0, 1]  
"""
pre y0  : θ=0, dθ/dt=1
""" 


for dt in dt_values:
    t = np.arange(0, t_max, dt)
    
    for omega_b in omega_b_values:
        sol_RK4 = RungeKutta4Order(pendulum, y0, t, args=(eta, omega, A, omega_b))
        sol_Euler = EulerMethod(pendulum, y0, t, args=(eta, omega, A, omega_b))
        

        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 1, 1)
        plt.plot(t, sol_RK4[:, 0], label=f'RK4, ωb={omega_b}')
        plt.plot(t, sol_Euler[:, 0], '--', label=f'Euler, ωb={omega_b}')
        plt.xlabel('Čas t')
        plt.ylabel('θ(t)')
        plt.title(f'Výchylka kyvadla, dt={dt}, ωb={omega_b}')
        plt.legend()
        plt.grid(True)
        
        plt.subplot(2, 1, 2)
        plt.plot(t, sol_RK4[:, 1], label=f'RK4, ωb={omega_b}')
        plt.plot(t, sol_Euler[:, 1], '--', label=f'Euler, ωb={omega_b}')
        plt.xlabel('Čas t')
        plt.ylabel('dθ/dt')
        plt.title(f'Úhlová rychlost, dt={dt}, ωb={omega_b}')
        plt.legend()
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig(f'pendulum_omega_b_{omega_b}_dt_{dt}.png')
        plt.show()

#Uloha 2 

def RungeKutta4Order(f, y0, t_span, dt, args=()):
    """Optimalizovaná RK4 metoda"""
    t_start, t_end = t_span
    n_steps = int((t_end - t_start) / dt)
    

    t_poincare = np.array([t_start + i*2*np.pi/args[3] for i in range(1, N+1)])
    poincare_indices = np.array([int(t/dt) for t in t_poincare])
    
    t = t_start
    y = np.array(y0, dtype=float)
    
    """Vytvorime nove pole na vnkladanie, namiesto ukladania celej trajektorie ukladame iba poincare body"""
    poincare_points = np.zeros((len(poincare_indices), len(y0)))
    poincare_idx = 0
    
  
    for i in range(n_steps):
        if poincare_idx < len(poincare_indices) and i == poincare_indices[poincare_idx]:
            poincare_points[poincare_idx] = y
            poincare_idx += 1
            
        k1 = f(y, t, *args)
        k2 = f(y + k1*dt/2, t + dt/2, *args)
        k3 = f(y + k2*dt/2, t + dt/2, *args)
        k4 = f(y + k3*dt, t + dt, *args)
        
        y = y + (k1 + 2*k2 + 2*k3 + k4)*dt/6
        t += dt
    
    return poincare_points


omega_b = 0.6  
T_b = 2*np.pi/omega_b  
dt = 0.01  
N = 1000 


initial_conditions = [[0, omega0] for omega0 in np.linspace(0, 1, 10)]

plt.figure(figsize=(10, 8))


colors = plt.cm.jet(np.linspace(0, 1, len(initial_conditions)))

for i, y0 in enumerate(initial_conditions):
    poincare_points = RungeKutta4Order(
        pendulum, y0, [0, (N+1)*T_b], dt, 
        args=(eta, omega, A, omega_b)
    )
    
  
    poincare_theta = poincare_points[:, 0]
    poincare_omega = poincare_points[:, 1]
    

    plt.scatter(poincare_theta, poincare_omega, s=15, color=colors[i],
                label=f'$\\dot{{\\theta}}_0$={y0[1]:.1f}')

plt.xlabel('$\\theta(nT_b)$')
plt.ylabel('$\\dot{\\theta}(nT_b)$')
plt.title(f'Poincarého řez pro $\\omega_b$={omega_b}, $\\eta$={eta}, A={A}')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig(f'poincare_section_omega_b_{omega_b}.png')
plt.show()

#Poznamka k ulohe: Ked som googlil poincareho rez, tak som zistil, ze system ktory nam vypluje tato medota by mal byt "slabo chaoticky" co znamena ze akoze system je deterministicky, teda
#taky je aj u samotneho chaosu ale neprekracuje bifurkacne body. Akoze po spusteni pre 1000 bodov to skutocne vyzera ze to neni chaoticke ale ten vypočet trval 6 min cize mam asi nieco neefiktivne v kode, cize mozno odporucam spustit pre 100 bodov
