import numpy as np
from matplotlib import pyplot as plt 
import astro_tools

# Constants 
R = 8314.4621/44.1
gamma = 1.37

# Mars standard atmosphere 

def get_pressure(h):

    if h < 39000:
        p = 610.5*(228.5/(228.5-1.8*h/1000))**(19.435/-1.8)
    elif h < 48000:
        p = 11.6025*np.exp(-19.435*(h/1000-39)/158.3)
    elif h < 55000:
        p = 3.84305*(158.3/(158.3-2.35*(h/1000-48)))**(19.435/-2.35)
    elif h < 66000:
        p = 1.55091*(141.85/(141.85+0.65*(h/1000-55)))**(19.435/0.65)
    elif h < 75000:
        p = 0.356464*(149/(149-2.5*(h/1000-66)))**(19.435/-2.5)
    elif h < 84000:
        p = 0.099843*(126.5/(126.5+2.5*(h/1000-75)))**(19.435/2.5)
    elif h < 95000:
        p = 0.0279653*np.exp(-19.435*(h/1000-84)/149)
    elif h < 105000:
        p = 0.00666032*(149/(149-1.4*(h/1000-95)))**(19.435/-1.4)
    elif h < 115897:
        p = 0.00169282*(135/(135-0.65*(h/1000-105)))**(19.435/-0.65)
    else:
        p = 0
    return p    

def get_temperature(h):

    if h < 39000:
        T = 228.5 - 1.8*h/1000
    elif h < 48000:
        T = 158.3
    elif h < 55000:
        T = 271.1 - 2.35*h/1000
    elif h < 66000:
        T = 106.1 + 0.65*h/1000
    elif h < 75000:
        T = 314 - 2.5*h/1000
    elif h < 84000:
        T = -61 + 2.5*h/1000
    elif h < 95000:
        T = 149
    elif h < 105000:
        T = 282 - 1.4*h/1000
    elif h < 115897:
        T = 203.25 - 0.65*h/1000
    else:
        e = (h/1000 - 120)*(3389.51 + 120)/(3389.51 + h/1000)
        T = 200 - 72.225*np.exp(-0.0195*e)
    return T   

def get_density(p, T, R=R):
    return p/(R*T)

def get_speed_of_sound(T, gamma=gamma, R=R):
    return np.sqrt(gamma*R*T)

if __name__=="__main__":
    alt = np.arange(0,100000,1)
    pressure = [get_pressure(h) for h in alt]
    temp = [get_temperature(h) for h in alt]
    density = [get_density(get_pressure(h),get_temperature(h)) for h in alt]
    mars = astro_tools.Planet()
    rho = [mars.density(h) for h in alt]
        
    plt.rcParams.update({'font.size': 12})

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Altitude [km]')
    ax1.set_ylabel('Temperature [K]', color=color)
    ax1.plot(np.array(alt)/1000, temp, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    plt.grid()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Pressure [Pa]', color=color)  # we already handled the x-label with ax1
    ax2.plot(np.array(alt)/1000, pressure, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    plt.plot(alt, density, label = 'Density from Viking I entry data')
    plt.plot(alt,rho, label = 'Exponentail density model')
    plt.grid()
    plt.legend(loc = 'best')
    plt.ylabel('Density [kg/m^3]')
    plt.xlabel('Altitude [km]')
    plt.show()