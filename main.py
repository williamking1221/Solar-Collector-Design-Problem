import numpy as np
import matplotlib.pyplot as plt


class SolarCollector:
    def __init__(self, collector_type):
        """
        Initializes Solar Collector, and fills in Properties from type given
        :param collector_type: Solar Collector Type: 1 == DF 100-10, 2 == DF 100-20, 3 == DF 100-30
        """
        self.collector_type = collector_type
        if self.collector_type == 1:
            self.aperture_area = 1.07  # m^2
        elif self.collector_type == 2:
            self.aperture_area = 2.15  # m^2
        elif self.collector_type == 3:
            self.aperture_area = 3.23  # m^2
        else:
            print("Error: Type must be either 1, 2, or 3")
        self.nu = 0.773  # nu_0
        self.alpha_1a = 1.43  # W/(m^2K)
        self.alpha_2a = 0.006  # W/(m^2K^2)
        self.flow_rate = 0.020 * self.aperture_area  # kg/s


class Container:
    def __init__(self, height, radius):
        """
        Initializes Container without insulation
        :param height: vertical height of container (m)
        :param radius: radius of container (m)
        """
        self.height = height
        self.radius = radius
        self.volume = (np.pi * radius ** 2) * height
        self.insulation_resistance = None
        self.insulation_thickness = None
        self.insulation_k_val = None
        self.k_val = None
        self.wall_temp = None
        self.outside_h = None
        self.u_val = None

    def set_insulation(self, k, thickness):
        """
        Set an insulation layer to the container, with a given k value
        :param k: Conduction Coefficient (W/mK)
        :param thickness: thickness of the insulation (m)
        :return: None
        """
        self.k_val = k
        self.insulation_thickness = thickness
        if thickness > 0:
            self.insulation_resistance = np.log((thickness + self.radius) / self.radius) / (2 * np.pi * self.height * k)

    def set_outside_temp(self, outside_temp, **kwargs):
        """
        Set wall temperature of the container, assumed to be 15 degrees Celsius higher than the outside temperature
        :param outside_temp: Ambient air temperature outside the container (K)
        :param kwargs: Properties of air at the outside temp, such as alpha, beta, viscosity
        :return: None
        """
        self.wall_temp = outside_temp + 15
        keys = kwargs.keys()
        if "beta" in keys:
            beta = kwargs["beta"]
        else:
            beta = 1 / outside_temp  # K^-1
        if "alpha" in keys:
            alpha = kwargs["alpha"]  # K^-1
        else:
            alpha = 0
            print("No Alpha, Default estimate of h = 10 W/m^2K applied")
            self.outside_h = 10  # W/m^2K
            return
        if "viscosity" in keys:
            viscosity = kwargs["viscosity"]  # N * s /m^2
        else:
            viscosity = 0
            print("No Viscosity, Default estimate of h = 10 W/m^2K applied")
            self.outside_h = 10  # W/m^2K
            return
        if self.k_val is None or self.k_val == 0:
            print("No k_val, Default estimate of h = 10 W/m^2K applied")
            self.outside_h = 10  # W/m^2K
            return
        g = 9.81  # m/s^2
        delta_t = 15  # K
        rayleigh_num = (g * beta * delta_t * self.height ** 3) / (viscosity * alpha)
        if rayleigh_num > 10 ** 9:
            self.outside_h = (self.k_val * 0.13 * rayleigh_num ** (1 / 3)) / self.height
            print("Outside h = {:.2f}".format(self.outside_h))
        else:
            if "prandtl" in keys:
                p_num = kwargs["prandtl"]
                nusselt_num = (p_num / (p_num + 0.986 * p_num ** (1 / 2) + 0.492)) ** (1 / 4) * rayleigh_num ** (1 / 4)
                self.outside_h = nusselt_num * self.k_val / self.height
                print("Outside h = {:.2f}".format(self.outside_h))
            else:
                p_num = 0
                print("No Prandtl, Default estimate of h = 10 W/m^2K applied")
                self.outside_h = 10  # W/m^2K
                return

    def optimize_insulation_thickness(self):
        if self.k_val is None or self.k_val == 0:
            print("No Optimization Required, no k_val given")
        if self.outside_h is None or self.outside_h == 0:
            print("No Optimization Possible, no h_value calculated yet")
        self.insulation_thickness = self.k_val / self.outside_h     # m
        print('New Optimized Insulation Thickness is {:.2f} cm'.format(self.insulation_thickness * 100))

    def calculate_u_val(self):
        if self.k_val is None and self.outside_h is None:
            print("Invalid input -- No values for outside h, or for k")
            return
        if self.k_val is None or self.k_val == 0:
            self.u_val = self.outside_h
            return
        else:
            self.u_val = ((1 / self.outside_h) * (self.radius / (self.radius + self.insulation_thickness)) +
                          (self.radius / self.k_val) * np.log((self.radius + self.insulation_thickness) / self.radius))
            return


water_cp_dict = {280: 4.198, 283: 4.198 + (3 / 5) * (4.189 - 4.198), 285: 4.189, 290: 4.184, 295: 4.181, 300: 4.179,
                 305: 4.178, 310: 4.178, 315: 4.179, 320: 4.180, 325: 4.182, 330: 4.184, 335: 4.186, 340: 4.188,
                 345: 4.191, 350: 4.195}  # kJ/kgK

water_density_dict = {280: 1, 283: 1, 285: 1, 290: 0.999, 295: 0.998, 300: 0.997, 305: 0.995, 310: 0.993, 315: 0.991,
                      320: 0.989, 325: 0.987, 330: 0.984, 335: 0.982, 340: 0.979, 345: 0.9765, 350: 0.9737}  # kg / L

water_consumption_dict = {6: 0.2, 7: 3.8, 8: 15.4, 9: 5.3, 10: 0.5, 11: 4.9, 12: 7.8, 13: 3.6, 14: 0.4, 15: 0.2,
                          16: 3.4, 17: 7.6, 18: 8.5, 19: 12.3, 20: 4.2, 21: 3.7, 22: 0.3, 23: 0}  # gal / hr

solar_irradiance_dict = {6: 12.1, 7: 44.5, 8: 197.2, 9: 462.5, 10: 550.3, 11: 667.6, 12: 688.6, 13: 771.5, 14: 585.5,
                         15: 516.0, 16: 371.3, 17: 71.8, 18: 56.4, 19: 27.7, 20: 2.9, 21: 0, 22: 0, 23: 0}  # W/m^2


def get_water_prop_estimate(temp, dic):
    """
    Linear Interpolation of a property of water, given a dictionary of that property
    :param temp: Temperature in question
    :param dic: Dictionary of property
    :return: Linearly interpolated value of dict for temperature of water
    """
    temp_keys = sorted(dic.keys())
    if temp < temp_keys[0] or temp > temp_keys[-1]:
        raise ValueError("Temperature is outside the range of provided data.")

    # Find the closest lower and upper temperature keys
    lower_temp = max(key for key in temp_keys if key <= temp)
    upper_temp = min(key for key in temp_keys if key >= temp)

    # Linear interpolation formula
    lower_value = dic[lower_temp]
    upper_value = dic[upper_temp]

    if upper_temp == lower_temp and upper_value == lower_value:
        return upper_value

    prop_estimate = lower_value + (upper_value - lower_value) * (temp - lower_temp) / (upper_temp - lower_temp)

    return prop_estimate


def solar_heat_solver(temperature_in, collector, irradiance):
    """
    A two iteration method of computing the roots of the quadratic equation generated by equating the two methods of
    computing for q_solar,gen
    :param temperature_in: Incoming temperature to the solar collector heat exchanger
    :param collector: collector, of class SolarCollector
    :param irradiance: irradiance at that hour, W/m^2
    :return: Heat Gain, Final temperature of the outflowing fluid
    """
    # First Iteration (Assume cp_temperature_in = cp used)
    cp_temperature_in = get_water_prop_estimate(temperature_in, water_cp_dict) * 1000  # J/kgK
    a = collector.alpha_2a * collector.aperture_area
    b_first_iter = collector.alpha_1a * collector.aperture_area + collector.flow_rate * cp_temperature_in
    c = - collector.aperture_area * collector.nu * irradiance
    deltas_first_iter = np.roots(np.array([a, b_first_iter, c]))  # K
    # Take real roots
    real_deltas_first_iter = [root.real for root in deltas_first_iter if np.isreal(root)]
    if real_deltas_first_iter:
        delta_first_iter = max(real_deltas_first_iter)
    else:
        print("Warning No Real Roots")
        delta_first_iter = 0
    # Second Iteration
    temperature_out = temperature_in + delta_first_iter
    cp_temperature_out = get_water_prop_estimate(temperature_out, water_cp_dict) * 1000
    cp_next_iter = (cp_temperature_out + cp_temperature_in) / 2  # Take average
    b_next_iter = collector.alpha_1a * collector.aperture_area + collector.flow_rate * cp_next_iter
    deltas_next_iter = np.roots(np.array([a, b_next_iter, c]))  # K
    # Take real roots
    real_deltas_next_iter = [root.real for root in deltas_next_iter if np.isreal(root)]
    if real_deltas_next_iter:
        delta_next_iter = max(real_deltas_next_iter)
    else:
        print("Warning No Real Roots")
        delta_next_iter = 0
    # Final Values
    temperature_out_final = temperature_in + delta_next_iter
    cp_out_final = get_water_prop_estimate(temperature_out_final, water_cp_dict) * 1000
    cp_average_final = (cp_out_final + cp_temperature_in) / 2
    total_q_solar_gen = delta_next_iter * collector.flow_rate * cp_average_final * 3600  # J
    return total_q_solar_gen, temperature_out_final


def home_heat_loss(temperature_in, flow_rate):
    """
    Determine heat loss for a certain hour
    :param temperature_in: Temperature of the water in the tank (K)
    :param flow_rate: Water Consumption for the hour (gal / hr)
    :return: heat loss of the process (J)
    """
    cp_temperature_in = get_water_prop_estimate(temperature_in, water_cp_dict) * 1000
    cp_temperature_out = water_cp_dict[283]
    cp_average = (cp_temperature_out + cp_temperature_in) / 2
    heat_loss = flow_rate * 3.785 * cp_average * (temperature_in - 283.5)  # J
    return heat_loss


def container_heat_loss(current_water_temp, outside_temp, container):
    """
    Calculate heat loss due to natural convection
    :param current_water_temp: Current Temperature of water in the container (K)
    :param outside_temp: Current Outside Temperature (K)
    :param container: Container, with u-value
    :return: Heat Loss (J)
    """
    surface_area = (np.pi * container.radius * 2) * container.height  # m^2
    temperature_diff = current_water_temp - outside_temp  # K
    heat_loss = container.u_val * temperature_diff * surface_area * 3600  # J
    return heat_loss


def calculate_new_water_temp(heat_loss_from_container, heat_loss_from_home, solar_heat_gain, current_water_temp,
                             current_water_volume):
    """
    Calculate the new water temperature in the tank
    :param heat_loss_from_container: Container heat loss (see method container_heat_loss above) (J)
    :param heat_loss_from_home: Home Water Consumption heat loss (see method home_heat_loss above) (J)
    :param solar_heat_gain: Heat Gain from Solar Collector (see method solar_heat_solver above) (J)
    :param current_water_temp: Temperature of Water in Container (K)
    :param current_water_volume: Volume of Water in Container (L)
    :return: New Water Temperature in Container (K)
    """
    net_heat_gain = solar_heat_gain - heat_loss_from_home - heat_loss_from_container  # J
    cp_estimate = get_water_prop_estimate(current_water_temp, water_cp_dict) * 1000  # J/kg*K
    density_estimate = get_water_prop_estimate(current_water_temp, water_density_dict)  # kg/L
    temp_diff = net_heat_gain / (cp_estimate * current_water_volume * density_estimate)  # K
    return current_water_temp + temp_diff


def simulate(tank_radius, tank_height, insulation_thickness, insulation_k, initial_water_temp, water_volume,
             num_type_1_solar_collector=0, num_type_2_solar_collector=0, num_type_3_solar_collector=0,
             outside_air_temp=300, optimize=True, **kwargs):
    """
    Simulate hour by hour heat loss and temperature change
    :param tank_radius: Radius of the Container (m)
    :param tank_height: Height of the Container (m)
    :param insulation_thickness: Thickness of Insulation for the container (m)
    :param insulation_k: k-value of the insulation for the container (W/mK)
    :param initial_water_temp: Initial Temperature of the Water at hour 0 (K)
    :param water_volume: Initial Volume of the Water in the tank (L)
    :param num_type_1_solar_collector: Number of Type 1 Solar Collector
    :param num_type_2_solar_collector: Number of Type 2 Solar Collector
    :param num_type_3_solar_collector: Number of Type 3 Solar Collector
    :param outside_air_temp: Outside air temperature (K), fixed at 300
    :param optimize: Optimize insulation radius, default True
    :param kwargs: alpha, beta, viscosity, and prandtl numbers to determine convection coefficient for u-value calculation
    :return: None, but prints statements on details
    """
    water_tank = Container(tank_height, tank_radius)
    water_tank.set_insulation(insulation_k, insulation_thickness)
    water_tank.set_outside_temp(outside_air_temp, **kwargs)
    if optimize:
        water_tank.optimize_insulation_thickness()
    water_tank.calculate_u_val()
    solar_collectors = []

    if (num_type_1_solar_collector + num_type_2_solar_collector + num_type_3_solar_collector) == 0:
        print("Not Feasible, add at least one solar collector")
        return
    for i in range(num_type_1_solar_collector):
        solar_collectors.append(SolarCollector(1))
    for i in range(num_type_2_solar_collector):
        solar_collectors.append(SolarCollector(2))
    for i in range(num_type_3_solar_collector):
        solar_collectors.append(SolarCollector(3))

    if water_volume > (tank_height * np.pi * (tank_radius ** 2) * 1000):
        print("Not Feasible, need a larger tank")
    minimum_water_level_req = 3.785 * 15.4 + 0.020 * (num_type_1_solar_collector * 1.07 * 3600 +
                                                      num_type_2_solar_collector * 2.15 * 3600 +
                                                      num_type_3_solar_collector * 3.23 * 3600)    # L
    if water_volume < minimum_water_level_req:
        print("Not Feasible, Will Lack Water at some time. Current Water Requirement is {:.2f} L"
              .format(minimum_water_level_req))
        return

    current_water_temp = initial_water_temp
    current_water_volume = water_volume
    final_temperatures = []
    solar_collector_heat_gain = []
    home_consumption_heat_loss = []
    container_natural_heat_loss = []
    net_heat_gain = []

    # Hour 0 - 5
    for i in range(6):
        container_heat_loss_for_hour = container_heat_loss(current_water_temp, outside_air_temp, water_tank)
        total_heat_gain_for_hour = -container_heat_loss_for_hour
        new_water_temp = calculate_new_water_temp(container_heat_loss_for_hour, 0, 0, current_water_temp,
                                                  current_water_volume)
        print("Hour: {} || Heat Gain from Solar Collectors: 0.00 J || Heat Loss from Home Consumption: 0.00 J || Heat "
              "Loss from Container: {:.2f} J || Net Heat Gain: {:.2f} J ||New Temperature: {:.2f} K"
              .format(i, container_heat_loss_for_hour, total_heat_gain_for_hour, new_water_temp))
        current_water_temp = new_water_temp
        final_temperatures.append(current_water_temp)
        solar_collector_heat_gain.append(0)
        home_consumption_heat_loss.append(0)
        container_natural_heat_loss.append(container_heat_loss_for_hour / 1000)
        net_heat_gain.append(total_heat_gain_for_hour / 1000)
    for i in solar_irradiance_dict.keys():
        solar_heat_gain_for_hour = 0
        home_heat_loss_for_hour = 0
        container_heat_loss_for_hour = 0
        for j in solar_collectors:
            solar_heat_gain_for_hour += solar_heat_solver(current_water_temp, j, solar_irradiance_dict[i])[0]
        home_heat_loss_for_hour = home_heat_loss(current_water_temp, water_consumption_dict[i])
        container_heat_loss_for_hour = container_heat_loss(current_water_temp, outside_air_temp, water_tank)
        total_heat_gain_for_hour = solar_heat_gain_for_hour - home_heat_loss_for_hour - container_heat_loss_for_hour
        new_water_temp = calculate_new_water_temp(container_heat_loss_for_hour, home_heat_loss_for_hour,
                                                  solar_heat_gain_for_hour, current_water_temp, current_water_volume)
        print("Hour: {} || Heat Gain from Solar Collectors: {:.2f} J || Heat Loss from Home Consumption: {:.2f} J || "
              "Heat Loss from Container: {:.2f} J || Net Heat Gain: {:.2f} J ||New Temperature: {:.2f} K"
              .format(i, solar_heat_gain_for_hour, home_heat_loss_for_hour, container_heat_loss_for_hour,
                      total_heat_gain_for_hour, new_water_temp))
        current_water_temp = new_water_temp
        final_temperatures.append(current_water_temp)
        solar_collector_heat_gain.append(solar_heat_gain_for_hour / 1000)
        home_consumption_heat_loss.append(home_heat_loss_for_hour/ 1000)
        container_natural_heat_loss.append(container_heat_loss_for_hour / 1000)
        net_heat_gain.append(total_heat_gain_for_hour / 1000)

    # Plotting
    hour_scale = np.linspace(0, 24, 24)
    fig1 = plt.figure(num=1, clear=True)
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.plot(hour_scale, final_temperatures)
    ax1.grid(True)
    ax1.set_xlabel("Hour")
    ax1.set_ylabel("Final Temperature at End of Hour (K)")
    ax1.set_title("Plot of Final Temperature vs Hour")
    fig1.savefig("Plot of Final Temperature vs Hour")

    fig2 = plt.figure(num=2, clear=True)
    ax2 = fig2.add_subplot(1, 1, 1)
    ax2.plot(hour_scale, solar_collector_heat_gain)
    ax2.grid(True)
    ax2.set_xlabel("Hour")
    ax2.set_ylabel("Heat Gain from Solar Collectors (kJ)")
    ax2.set_title("Plot of Final Solar Collector Heat Gain vs Hour")
    fig2.savefig("Plot of Final Solar Collector Heat Gain vs Hour")

    fig3 = plt.figure(num=3, clear=True)
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.plot(hour_scale, home_consumption_heat_loss)
    ax3.grid(True)
    ax3.set_xlabel("Hour")
    ax3.set_ylabel("Heat Loss from House Consumption (kJ)")
    ax3.set_title("Plot of House Consumption Heat Loss vs Hour")
    fig3.savefig("Plot of House Consumption Heat Loss vs Hour")

    fig4 = plt.figure(num=4, clear=True)
    ax4 = fig4.add_subplot(1, 1, 1)
    ax4.plot(hour_scale, container_natural_heat_loss)
    ax4.grid(True)
    ax4.set_xlabel("Hour")
    ax4.set_ylabel("Heat Loss within Container due to Outside Natural Convection (kJ)")
    ax4.set_title("Plot of Container Heat Loss vs Hour")
    fig4.savefig("Plot of Container Heat Loss vs Hour")

    fig5 = plt.figure(num=5, clear=True)
    ax5 = fig5.add_subplot(1, 1, 1)
    ax5.plot(hour_scale, net_heat_gain)
    ax5.grid(True)
    ax5.set_xlabel("Hour")
    ax5.set_ylabel("Net Heat Gain in the Container (kJ)")
    ax5.set_title("Net Heat Gain vs Hour")
    fig5.savefig("Net Heat Gain vs Hour")
    plt.show()
    return


if __name__ == "__main__":
    container_tank_radius = 1       # m
    container_tank_height = 1       # m
    test_insulation_thickness = 0.077        # m (Optimized Thickness is 0.063 m)
    test_insulation_k = 0.046            # W/mK, for Mineral wool granules Foam
    test_initial_water_temp = 330        # K
    test_water_volume = (container_tank_height * np.pi * container_tank_radius ** 2) * 1000                # L
    type_1_solar_collectors = 1
    type_2_solar_collectors = 0
    type_3_solar_collectors = 1
    test_outside_air_temp = 300          # K
    beta_val = 1 / test_outside_air_temp   # K^-1
    alpha_val = 22.5 * (10 ** -6)   # m^2/s
    prandtl_val = 0.707
    viscosity_val = 15.89 * (10 ** -6)  # m^2/s

    simulate(container_tank_radius, container_tank_height, test_insulation_thickness, test_insulation_k,
             test_initial_water_temp, test_water_volume, num_type_1_solar_collector=type_1_solar_collectors,
             num_type_2_solar_collector=type_2_solar_collectors, num_type_3_solar_collector=type_3_solar_collectors,
             outside_air_temp=test_outside_air_temp, optimize=False,
             beta=beta_val, alpha=alpha_val, prandtl=prandtl_val, viscosity=viscosity_val)
