import numpy as np
import matplotlib.pyplot as plt
import os
import json as js


class MechaikCalculator:

    def __init__(self, calculations_steps) -> None:
        
        self.log_dir_path = "C:\\Users\\1\\Desktop\\MSF_lab_02\\log_dir"
        self.indignant_log_path = os.path.join(self.log_dir_path, "indignant_params.json")
        self.none_indignant_log_path = os.path.join(self.log_dir_path, "none_indignant_params.json")

        self.calculations_steps = calculations_steps
        self.virtual_time = 0.01

        self.H_apogee = 740 + 6371
        self.H_peregee = 350 + 6371

        self.semimajor_axis = (self.H_apogee + self.H_peregee) / 2.0
        self.ext_param = (self.H_apogee + self.H_peregee) / (self.H_apogee - self.H_peregee)
        self.focal_param = (self.H_apogee - self.H_peregee) / (2.0 * self.semimajor_axis)

        self.Omega_start = (40.0 * np.pi) / 180.0
        self.u_start = (70.0 * np.pi) / 180.0
        self.i_start = (37.1 * np.pi) / 180.0
        self.Sa_param = 23.0
        self.fly_mass = 2700
        self.c_coeff = 3.5

        self.earth_grav_param = 398600.4415
        self.atmosphere_density = 2.77744

        self.curent_radius_param = 0.0
        self.curent_u_param = self.u_start
        self.curent_Omega_param = self.Omega_start
        self.curent_omega_param = 0
        self.curent_i_param = self.i_start
        self.curent_focal_param = self.focal_param
        self.curent_ext_param = self.ext_param
    

    def _calculate_params(self, theta, simulation_type):

        self.curent_radius_param = self.curent_focal_param / (1.0 + self.curent_ext_param * np.cos(theta))
        
        if simulation_type == "indignant":
        

            self.velocity = (self.earth_grav_param / self.curent_focal_param) * (1 + (self.curent_ext_param ** 2) + (2 * self.curent_ext_param * np.cos(theta)))
            self.force_A = (self.c_coeff * self.atmosphere_density * self.velocity * self.Sa_param) / (2.0 * self.fly_mass)
            self.force_1 = -self.force_A * ((self.curent_ext_param * np.sin(theta)) / (np.sqrt(1 + (self.curent_ext_param ** 2) + (2.0 * self.curent_ext_param * np.cos(theta)))))
            self.force_2 = -self.force_A * ((1.0 + self.curent_ext_param * np.cos(theta)) / (np.sqrt(1 + (self.curent_ext_param ** 2) + (2.0 * self.curent_ext_param * np.cos(theta)))))
        
        elif simulation_type == "none_indignant":

            self.force_1 = 0
            self.force_2 = 0

        delta_focal_param = ((2 * self.curent_radius_param ** 2) / self.earth_grav_param) * self.force_2
        delta_ext_param = ((self.curent_radius_param ** 2) / self.earth_grav_param) * (self.force_1 * np.sin(theta) + self.force_2 * ((1 + self.curent_radius_param / self.curent_focal_param) * np.cos(theta) + self.curent_ext_param * (self.curent_radius_param / self.curent_focal_param)))
        delta_omega_param = ((self.curent_radius_param ** 2) / (self.earth_grav_param * self.curent_ext_param)) * (-self.force_1 * np.cos(theta) + self.force_2 * (1 + (self.curent_radius_param / self.curent_focal_param) * np.sin(theta)))

        delta_focal_param = (delta_focal_param * np.pi) / 180.0
        delta_ext_param = (delta_ext_param * np.pi) / 180.0
        delta_omega_param = (delta_omega_param * np.pi) / 180.0
        
        self.theta_data.append(self.hiden_theta)
        self.radius_data.append(self.curent_radius_param)
        self.Omega_data.append(self.curent_Omega_param)
        self.omega_data.append(self.curent_omega_param)
        self.focal_data.append(self.curent_focal_param)
        self.ext_data.append(self.curent_ext_param)
        self.i_data.append(self.curent_i_param)
        self.u_data.append(self.curent_u_param)

        self.curent_ext_param += delta_ext_param * self.virtual_time
        self.curent_focal_param += delta_focal_param * self.virtual_time
        self.curent_omega_param += delta_omega_param * self.virtual_time
        self.curent_Omega_param += 0
        self.curent_i_param += 0
        self.curent_u_param += 0

    

    def _save_data(self, dtype):

        data_json_format = {}
        print(len(self.Omega_data))
        print(len(self.omega_data))
        print(len(self.i_data))
        print(len(self.u_data))
        print(len(self.focal_data))
        print(len(self.ext_data))
        print(len(self.radius_data))
        print(len(self.theta_data))
        for iteration in range(self.calculations_steps):

            data_json_format[f"iteration_number: {iteration}"] = {

                "Omega_param": self.Omega_data[iteration],
                "omega_param": self.omega_data[iteration],
                "i_param": self.i_data[iteration],
                "u_param": self.u_data[iteration],
                "ext_param": self.ext_data[iteration],
                "focal_param": self.focal_data[iteration],
                "radius_param": self.radius_data[iteration],
                "theta_param": self.theta_data[iteration]
            }

        if dtype == "indignant":

            with open(self.indignant_log_path, "w") as file:
                js.dump(data_json_format, file)
            
        else:

            with open(self.none_indignant_log_path, "w") as file:
                js.dump(data_json_format, file)

    def show_data(self):

        omega_data_ni = []
        Omega_data_ni = []
        u_data_ni = []
        i_data_ni = []
        focal_data_ni = []
        ext_data_ni = []
        rad_data_ni = []
        theta_data_ni = []

        omega_data_i = []
        Omega_data_i = []
        u_data_i = []
        i_data_i = []
        focal_data_i = []
        ext_data_i = []
        rad_data_i = []
        theta_data_i = []
        


        with open(self.indignant_log_path, "r") as json_file:

            data_json = js.load(json_file)
            for iteration_step_key in data_json.keys():
                
                Omega_data_i.append(data_json[iteration_step_key]["Omega_param"])
                omega_data_i.append(data_json[iteration_step_key]["omega_param"])
                u_data_i.append(data_json[iteration_step_key]["u_param"])
                i_data_i.append(data_json[iteration_step_key]["i_param"])
                focal_data_i.append(data_json[iteration_step_key]["focal_param"])
                ext_data_i.append(data_json[iteration_step_key]["ext_param"])
                rad_data_i.append(data_json[iteration_step_key]["radius_param"])
                theta_data_i.append(data_json[iteration_step_key]["theta_param"])

        with open(self.none_indignant_log_path, "r") as json_file:

            data_json = js.load(json_file)
            for iteration_step_key in data_json.keys():
                
                Omega_data_ni.append(data_json[iteration_step_key]["Omega_param"])
                omega_data_ni.append(data_json[iteration_step_key]["omega_param"])
                u_data_ni.append(data_json[iteration_step_key]["u_param"])
                i_data_ni.append(data_json[iteration_step_key]["i_param"])
                focal_data_ni.append(data_json[iteration_step_key]["focal_param"])
                ext_data_ni.append(data_json[iteration_step_key]["ext_param"])
                rad_data_ni.append(data_json[iteration_step_key]["radius_param"])
                theta_data_ni.append(data_json[iteration_step_key]["theta_param"])
                
        plt.style.use("dark_background")
        fig, axis = plt.subplots(nrows=2, ncols=3)

        axis[0, 0].plot(Omega_data_ni, color="g", label="Omega")
        axis[0, 1].plot(omega_data_ni, color="g", label="omega")
        axis[0, 2].plot(u_data_ni, color="g", label="u")
        axis[1, 0].plot(i_data_ni, color="g", label="i")
        axis[1, 1].plot(focal_data_ni, color="g", label="focal")
        axis[1, 2].plot(ext_data_ni, color="g", label="ext")

        axis[0, 0].legend(loc="upper left")
        axis[0, 1].legend(loc="upper left")
        axis[0, 2].legend(loc="upper left")
        axis[1, 0].legend(loc="upper left")
        axis[1, 1].legend(loc="upper left")
        axis[1, 2].legend(loc="upper left")

        plt.show()

        fig, axis = plt.subplots(nrows=2, ncols=3)

        axis[0, 0].plot(Omega_data_i, color="r", label="Omega")
        axis[0, 1].plot(omega_data_i, color="r", label="omega")
        axis[0, 2].plot(u_data_i, color="r", label="u")
        axis[1, 0].plot(i_data_i, color="r", label="i")
        axis[1, 1].plot(focal_data_i, color="r", label="focal")
        axis[1, 2].plot(ext_data_i, color="r", label="ext")

        axis[0, 0].legend(loc="upper left")
        axis[0, 1].legend(loc="upper left")
        axis[0, 2].legend(loc="upper left")
        axis[1, 0].legend(loc="upper left")
        axis[1, 1].legend(loc="upper left")
        axis[1, 2].legend(loc="upper left")
        
        plt.show()
        
        plt.style.use("dark_background")
        fig, axis = plt.subplots(subplot_kw={'projection': 'polar'}, ncols=2)
    
        axis[0].plot(theta_data_ni, rad_data_ni, color="g", label="none indignants trajectory")
        axis[1].plot(theta_data_i, rad_data_ni, color="r", label="indignant trajectory")
        axis[0].legend(loc="upper left")
        axis[1].legend(loc="upper left")

        plt.show()

    def run_simulation(self, simulation_mode="none_indignant"):
        

        self.u_data = []
        self.Omega_data = []
        self.omega_data = []
        self.i_data = []
        self.focal_data = []
        self.ext_data = []
        self.radius_data = []
        self.theta_data = []

        self.curent_radius_param = 0.0
        self.curent_u_param = self.u_start
        self.curent_Omega_param = self.Omega_start
        self.curent_omega_param = 0
        self.curent_i_param = self.i_start
        self.curent_focal_param = self.focal_param
        self.curent_ext_param = self.ext_param

        for theta in np.linspace(0, 3600 * np.pi, self.calculations_steps):
            
            self.hiden_theta = (theta * np.pi) / 180.0
            theta = (theta * np.pi) / 180.0
            self._calculate_params(simulation_type=simulation_mode, theta=theta)

        self._save_data(dtype=simulation_mode)



if __name__ == "__main__":

    sim_object = MechaikCalculator(calculations_steps=100000)
    sim_object.run_simulation()
    sim_object.run_simulation(simulation_mode="indignant")
    sim_object.show_data()


            

        
        


        
        