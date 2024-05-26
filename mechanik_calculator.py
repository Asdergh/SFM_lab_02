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
        self.virtual_time = 0.1

        self.H_apogee = 740 + 6371
        self.H_peregee = 350 + 6371

        self.semimajor_axis = (self.H_apogee + self.H_peregee) / 2.0
        self.ext_param = (self.H_apogee + self.H_peregee) / (self.H_apogee - self.H_peregee)
        self.focal_param = (self.H_apogee - self.H_peregee) / (2.0 * self.semimajor_axis)
        print(self.ext_param)
        print(self.semimajor_axis)
        print(self.focal_param)

        self.Gamma_start = (40.0 * np.pi) / 180.0
        self.u_start = (70.0 * np.pi) / 180.0
        self.i_start = (37.1 * np.pi) / 180.0
        self.Sa_param = 23.0
        self.fly_mass = 2700
        self.c_coeff = 3.5

        self.earth_grav_param = 398600.4415
        self.atmosphere_density = 2.77744

        theta = 0
        self.curent_u_param = 0
        self.curent_Omega_param = 0
        self.curent_omega_param = 0
        self.curent_i_param = 0
        self.curent_focal_param = 0
        self.curent_ext_param = 0
    

    def _calculate_none_indignant_SP(self, theta):

        self.radius = self.focal_param / (1.0 - self.ext_param * np.cos(theta))
        print(self.radius)
        self.gamma_inv = 1 / (1 - ((self.force_3 * self.radius ** 3) / (self.earth_grav_param * self.focal_param)) * np.sin(self.u_start) * (np.sin(self.u_start) / np.cos(self.u_start)))
        self.curent_u_param = self.u_start + (np.sqrt(self.earth_grav_param * self.focal_param) / self.radius ** 2) * self.gamma_inv * self.virtual_time
        self.curent_i_param = self.i_start + self.virtual_time * self.gamma_inv * (self.radius ** 3 / self.earth_grav_param * self.focal_param) * self.force_3 * np.cos(self.u_start)
        self.curent_focal_param = self.focal_param + self.virtual_time * self.gamma_inv * 2.0 * self.force_2 * (self.radius ** 3 / self.earth_grav_param)
        self.curent_Omega_param += self.virtual_time * self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.focal_param)) * self.force_3 * (np.sin(self.u_start) / np.sin(self.i_start)) 
        self.curent_omega_param += self.virtual_time * (self.radius ** 2 / (self.earth_grav_param * self.ext_param)) * (-self.force_1 * np.cos(theta) + self.force_2 * (1 + self.radius / self.focal_param) * np.sin(theta) - self.force_3 * (self.radius / self.focal_param) * self.ext_param * (np.cos(self.i_start) / np.sin(self.i_start)) * np.sin(self.u_start))
        self.curent_ext_param += self.virtual_time * self.gamma_inv * (self.radius ** 2 / self.earth_grav_param) * (self.force_1 * np.sin(theta) + self.force_2 * ((1 + self.radius / self.focal_param) * np.cos(theta) + self.ext_param * (self.radius / self.focal_param)))

    def _calculate_none_indignant_CP(self, theta):

        self.radius = self.curent_focal_param * (1 - self.ext_param * np.cos(theta))
        print(self.radius)
        self.curent_i_param += (((self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.curent_focal_param)) * self.force_3 * np.cos(self.curent_u_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_u_param += ((((np.sqrt(self.earth_grav_param * self.curent_focal_param) / self.radius ** 2) * self.gamma_inv) * self.virtual_time) * np.pi) / 180.0
        self.gamma_inv = 1.0 / (1.0 - ((self.force_3 * self.radius) / (self.earth_grav_param * self.curent_focal_param)) * np.sin(self.curent_u_param) * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)))
        self.curent_Omega_param += (((self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.curent_focal_param)) * self.force_3 * np.cos(self.curent_u_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_focal_param += (((self.gamma_inv * 2 * self.force_2 * (self.radius ** 3 / self.earth_grav_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_omega_param += (((self.gamma_inv * (self.radius ** 2 / (self.earth_grav_param * self.ext_param)) * (-self.force_1 * np.cos(theta) + self.force_2 * (1 + (self.radius / self.curent_focal_param)) * np.sin(theta) + self.force_3 * (self.radius / self.curent_focal_param) * self.curent_ext_param * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)) * np.sin(self.curent_i_param))) * self.virtual_time) * np.pi) / 180.0
        self.curent_ext_param += (((self.gamma_inv * (self.radius ** 2 / self.earth_grav_param) * (self.force_1 * np.sin(theta) + self.force_2 * ((1 + (self.radius / self.curent_focal_param)) * np.cos(theta) + self.curent_ext_param * (self.radius / self.curent_focal_param)))) * self.virtual_time) * np.pi) / 180.0

    def _calculate_indignant_SP(self, theta):

        self.radius = self.focal_param * (1 - self.ext_param * np.cos(theta))
        print(self.radius)
        self.vellocity = (self.earth_grav_param / self.focal_param) * (1.0 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta))
        self.q_coeff = (self.atmosphere_density * self.vellocity) / 2
        self.force_A = (self.c_coeff * self.q_coeff * self.Sa_param) / self.fly_mass
                    
        self.force_1 = -self.force_A * (self.ext_param * np.sin(theta) / np.sqrt(1 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta)))
        self.force_2 = -self.force_A * ((1 + self.ext_param * np.cos(theta)) / np.sqrt(1 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta)))
        self.force_3 = 0
                    
        self.curent_i_param += ((self.i_start * self.virtual_time) * np.pi) / 180.0
        self.curent_u_param += ((self.u_start * self.virtual_time) * np.pi) / 180.0
        self.gamma_inv = 1.0 / (1.0 - ((self.force_3 * self.radius) / (self.earth_grav_param * self.focal_param)) * np.sin(self.curent_u_param) * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)))
        self.curent_Omega_param += (((self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.focal_param)) * self.force_3 * np.cos(self.curent_u_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_focal_param += ((self.focal_param * self.virtual_time) * np.pi) / 180.0
        self.curent_omega_param += (((self.gamma_inv * (self.radius ** 2 / (self.earth_grav_param * self.ext_param)) * (-self.force_1 * np.cos(theta) + self.force_2 * (1 + (self.radius / self.focal_param)) * np.sin(theta) + self.force_3 * (self.radius / self.focal_param) * self.ext_param * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)) * np.sin(self.curent_i_param))) * self.virtual_time) * np.pi) / 180.0
        self.curent_ext_param += ((self.ext_param * self.virtual_time) * np.pi) / 180.0
    
    def _calculate_indignant_CP(self, theta):

        self.radius = self.curent_focal_param * (1 - self.ext_param * np.cos(theta))
        print(self.radius)
        self.vellocity = (self.earth_grav_param / self.focal_param) * (1.0 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta))
        self.q_coeff = (self.atmosphere_density * self.vellocity) / 2
        self.force_A = (self.c_coeff * self.q_coeff * self.Sa_param) / self.fly_mass

        self.force_1 = -self.force_A * (self.ext_param * np.sin(theta) / np.sqrt(1 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta)))
        self.force_2 = -self.force_A * ((1 + self.ext_param * np.cos(theta)) / np.sqrt(1 + self.ext_param ** 2 + 2 * self.ext_param * np.cos(theta)))
        self.force_3 = 0

        theta += ((((np.sqrt(self.earth_grav_param * self.curent_focal_param) / self.radius ** 2) + (1 / self.curent_ext_param) * np.sqrt(self.curent_focal_param / self.earth_grav_param) * (self.force_1 * np.cos(theta))) * self.virtual_time) * np.pi) / 180.0
        self.curent_i_param += (((self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.curent_focal_param)) * self.force_3 * np.cos(self.curent_u_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_u_param += ((((np.sqrt(self.earth_grav_param * self.curent_focal_param) / self.radius ** 2) * self.gamma_inv) * self.virtual_time) * np.pi) / 180.0
        self.gamma_inv = 1.0 / (1.0 - ((self.force_3 * self.radius) / (self.earth_grav_param * self.curent_focal_param)) * np.sin(self.curent_u_param) * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)))
        self.curent_Omega_param += (((self.gamma_inv * (self.radius ** 3 / (self.earth_grav_param * self.curent_focal_param)) * self.force_3 * np.cos(self.curent_u_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_focal_param += (((self.gamma_inv * 2 * self.force_2 * (self.radius ** 3 / self.earth_grav_param)) * self.virtual_time) * np.pi) / 180.0
        self.curent_omega_param += (((self.gamma_inv * (self.radius ** 2 / (self.earth_grav_param * self.ext_param)) * (-self.force_1 * np.cos(theta) + self.force_2 * (1 + (self.radius / self.curent_focal_param)) * np.sin(theta) + self.force_3 * (self.radius / self.curent_focal_param) * self.curent_ext_param * (np.cos(self.curent_i_param) / np.sin(self.curent_i_param)) * np.sin(self.curent_i_param))) * self.virtual_time) * np.pi) / 180.0
        self.curent_ext_param += (((self.gamma_inv * (self.radius ** 2 / self.earth_grav_param) * (self.force_1 * np.sin(theta) + self.force_2 * ((1 + (self.radius / self.curent_focal_param)) * np.cos(theta) + self.curent_ext_param * (self.radius / self.curent_focal_param)))) * self.virtual_time) * np.pi) / 180.0


    def _save_data(self, dtype):


        data_json_format = {}
        for (iteration_step, (Omega, omega, i, u, focal, ext, rad, theta)) in enumerate(zip(self.Omega_data_tensor, 
                                                                                self.omega_data_tensor,
                                                                                self.i_data_tensor,
                                                                                self.u_data_tensor,
                                                                                self.focal_data_tensor,
                                                                                self.ext_data_tensor,
                                                                                self.radius_data_tensor,
                                                                                self.theta_data_tensor)):
                
            data_json_format[f"iteration_step_number: {iteration_step}"] = {
                "Omega_param": Omega,
                "omega_param": omega,
                "i_param": i,
                "u_param": u,
                "focal_param": focal,
                "ext_param": ext,
                "rad_param": rad,
                "theta_param": theta,
            }

        if dtype == "indignant":

            self.curent_Omega_param = 0
            self.curent_omega_param = 0
            self.curent_i_param = 0
            self.curent_u_param = 0
            self.curent_focal_param = 0
            self.curent_ext_param = 0
            
            self.force_A = 0
            self.force_1 = 0
            self.force_2 = 0
            self.force_3 = 0
            self.q_coeff = 0

            with open(self.indignant_log_path, "w") as file:
                js.dump(data_json_format, file)


        elif dtype == "none indignant":
            
            self.curent_Omega_param = 0
            self.curent_omega_param = 0
            self.curent_i_param = 0
            self.curent_u_param = 0
            self.curent_focal_param = 0
            self.curent_ext_param = 0
            
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
                rad_data_i.append(data_json[iteration_step_key]["rad_param"])
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
                rad_data_ni.append(data_json[iteration_step_key]["rad_param"])
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

        print(theta_data_i, theta_data_ni, rad_data_i, rad_data_ni)
        axis[0].plot(rad_data_ni, theta_data_ni)
        axis[1].plot(rad_data_i, theta_data_i)

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

        if (simulation_mode == "none_indignant"):

            self.force_1 = 0
            self.force_2 = 0
            self.force_3 = 0
            
            for theta in np.arange(0, 10 * np.pi, 0.5):
                
                print(theta)
                if (theta == 0):
                    self._calculate_none_indignant_SP(theta)
                    
                else:
                    self._calculate_none_indignant_CP(theta)
                        
                self.Omega_data.append(self.curent_Omega_param)
                self.omega_data.append(self.curent_omega_param)
                self.i_data.append(self.curent_i_param)
                self.u_data.append(self.curent_u_param)
                self.focal_data.append(self.curent_focal_param)
                self.ext_data.append(self.curent_ext_param)
                self.radius_data.append(self.radius)
                self.theta_data.append(theta)

        
        elif (simulation_mode == "indignant"):

            for theta in np.arange(0, 10 * np.pi, 0.5):
                
                print(theta)
                if (theta == 0):
                    self._calculate_indignant_SP(theta)
                    

                else:
                    self._calculate_indignant_CP(theta)
                
                self.Omega_data.append(self.curent_Omega_param)
                self.omega_data.append(self.curent_omega_param)
                self.i_data.append(self.curent_i_param)
                self.u_data.append(self.curent_u_param)
                self.focal_data.append(self.curent_focal_param)
                self.ext_data.append(self.curent_ext_param)
                self.radius_data.append(self.radius)
                self.theta_data.append(theta)

        
        self.Omega_data_tensor = np.asarray(self.Omega_data)
        self.omega_data_tensor = np.asarray(self.omega_data)
        self.i_data_tensor = np.asarray(self.i_data)
        self.u_data_tensor = np.asarray(self.u_data)
        self.focal_data_tensor = np.asarray(self.focal_data)
        self.ext_data_tensor = np.asarray(self.ext_data)
        self.radius_data_tensor = np.asarray(self.radius_data)
        self.theta_data_tensor = np.asarray(self.theta_data)

        if simulation_mode == "indignant":
            self._save_data(dtype="indignant")
        
        elif simulation_mode == "none_indignant":
            self._save_data(dtype="none indignant")



if __name__ == "__main__":

    sim_object = MechaikCalculator(calculations_steps=100)
    sim_object.run_simulation()
    sim_object.run_simulation(simulation_mode="indignant")
    sim_object.show_data()


            

        
        


        
        