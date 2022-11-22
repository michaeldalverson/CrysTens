from typing import List
import os
import numpy as np
import math
from pymatgen.io.cif import CifParser

class CrysTens:
  def __init__(self, *args):
    if len(args) == 0:
      raise ValueError("Please input a path to a CIF file, a PyMatGen CIF Object, or the necessary parameters to create a CIF.")
    if len(args) == 9:
      self.from_params(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])
    elif isinstance(args[0], str):
      self.from_file(args[0])
    else:
      self.from_CIF_obj()
    
    self.direction_fix = 0.9999933
    self.normalize_near_max = {'atom': 83, 
                        'x': 0.979, 
                        'y': 0.9813, 
                        'z': 0.9782, 
                        'a': 19.866, 
                        'b': 26.43, 
                        'c': 42.5, 
                        'alpha': 107.69, 
                        'beta': 131.30000000000004, 
                        'gamma': 119.99999999999999, 
                        'sg': 227.0, 
                        'dir': 0.89146,
                        'length': 1.2779113271271993}

    self.normalize_avg_std = {'atom': [22.417031299395948, 27.503517675002012, 0.778701530199605, 2.7844790475732646],
                        'x': [0.2934632568241515, 0.42520262764739036, 0.6901727264665755, 1.661631491833322],
                        'y': [0.2919067377134818, 0.43638668432502087, 0.6689176095393177, 1.6226127447382617],
                        'z': [0.29308350406375905, 0.44429446323553035, 0.659660491668988, 1.591076536917126],
                        'a': [2.7537991622182245, 6.631296113153595, 0.14950307531158527, 6.018461573238658],
                        'b': [3.30605670472346, 6.887269725378788, 0.22413187899920117, 6.0803402458545905],
                        'c': [5.797509184209129, 8.88626989936126, 0.4084401245195147, 12.423378095203498],
                        'alpha': [2.7943744442379326, 89.96235708184788, 0, 1.3028296429486415],
                        'beta': [8.426823255100807, 92.73923379010695, 0, 1.4017387402522046],
                        'gamma': [12.722185599629745, 96.46287050653594, 0, 1.1121150950313212],
                        'sg': [76.31186396273776, 116.00356506238859, 0.6492202538967989, 1.3248570072359955],
                        'dir': [0.4099350434569409, 0.9999933333333323, 0.4099377763754505, 1.5900622236245516],
                        'length': [0.26288351681625777, 0.646245505143038, 0.4067858340586399, 2.15365474510736]}

    self.attribute_list = ["atom", "x", "y", "z", "a", "b", "c", "alpha", "beta", "gamma", "sg"]
    self.parameter_list = ["a", "b", "c", "alpha", "beta", "gamma", "sg"]



  def from_file(self, path: str):
    self.normalized = False
    parser = CifParser(path)
    structure = parser.get_structures()[0]
    self.a = structure.lattice.a
    self.b = structure.lattice.b
    self.c = structure.lattice.c

    self.alpha = structure.lattice.alpha
    self.beta = structure.lattice.beta
    self.gamma = structure.lattice.gamma

    self.sg = structure.get_space_group_info()[1]

    self.coord_list = []
    self.element_list = []
    for site in structure.sites:
      self.coord_list.append([site.a, site.b, site.c])
      self.element_list.append(site.specie.Z)
    
    self.construct_crys_tens()

  def from_params(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float, sg: int, element_list: List[str], coord_list :List[List[float]]):
    self.normalized = False
    self.a = a
    self.b = b
    self.c = c

    self.alpha = alpha
    self.beta = beta
    self.gamma = gamma

    self.sg = sg

    self.coord_list = coord_list
    self.element_list = element_list

    self.construct_crys_tens()

  def from_CIF_obj(self):
    #TODO
    pass

  def normalize_crys_tens(self, method: str = "divide_near_max"):
    if self.normalized:
      return
    self.normalized = True
    if method == "divide_near_max":
      for site_idx in range(self.crys_tens.shape[0] - 12):

        #Check if an atom is present
        if self.crys_tens[12 + site_idx, 0, 0] != 0.0:
          for att_idx, att in enumerate(self.attribute_list):
            self.crys_tens[att_idx, 12 + site_idx, :] /= self.normalize_near_max[att]
            self.crys_tens[12 + site_idx, att_idx, :] /= self.normalize_near_max[att]     

      for adj_idx in range(len(self.coord_list)):
        if site_idx != adj_idx:
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] /= self.normalize_near_max["length"]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] /= self.normalize_near_max["length"]

          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] += self.direction_fix
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] /= self.normalize_near_max["dir"]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] += self.direction_fix
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] /= self.normalize_near_max["dir"]

    elif method == "avg_std":
      for site_idx in range(self.crys_tens.shape[0] - 12):

        #Check if an atom is present
        if self.crys_tens[12 + site_idx, 0, 0] != 0.0:
          for att_idx, att in enumerate(self.attribute_list):
            self.crys_tens[att_idx, 12 + site_idx, :] -= self.normalize_avg_std[att][0]
            self.crys_tens[att_idx, 12 + site_idx, :] /= self.normalize_avg_std[att][1]
            self.crys_tens[att_idx, 12 + site_idx, :] += self.normalize_avg_std[att][2]
            self.crys_tens[att_idx, 12 + site_idx, :] /= self.normalize_avg_std[att][2] + self.normalize_avg_std[att][3]
            self.crys_tens[12 + site_idx, att_idx, :] -= self.normalize_avg_std[att][0]
            self.crys_tens[12 + site_idx, att_idx, :] /= self.normalize_avg_std[att][1]
            self.crys_tens[12 + site_idx, att_idx, :] += self.normalize_avg_std[att][2]
            self.crys_tens[12 + site_idx, att_idx, :] /= self.normalize_avg_std[att][2] + self.normalize_avg_std[att][3]   

      for adj_idx in range(len(self.coord_list)):

        if site_idx != adj_idx:
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] -= self.normalize_avg_std["length"][0]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] /= self.normalize_avg_std["length"][1]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] += self.normalize_avg_std["length"][2]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] /= self.normalize_avg_std["length"][2] + self.normalize_avg_std["length"][3]

          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] -= self.normalize_avg_std["length"][0]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] /= self.normalize_avg_std["length"][1]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] += self.normalize_avg_std["length"][2]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] /= self.normalize_avg_std["length"][2] + self.normalize_avg_std["length"][3]

          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] += self.direction_fix
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] -= self.normalize_avg_std["dir"][0]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] /= self.normalize_avg_std["dir"][1]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] += self.normalize_avg_std["dir"][2]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] /= self.normalize_avg_std["dir"][2] + self.normalize_avg_std["dir"][3]

          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] += self.direction_fix
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] -= self.normalize_avg_std["dir"][0]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] /= self.normalize_avg_std["dir"][1]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] += self.normalize_avg_std["dir"][2]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] /= self.normalize_avg_std["dir"][2] + self.normalize_avg_std["dir"][3]

  def unnormalize_crys_tens(self, method: str = "divide_near_max"):
    if not self.normalized:
      return
    self.normalized = False
    if method == "divide_near_max":
      for site_idx in range(self.crys_tens.shape[0] - 12):

        #Check if an atom is present
        if self.crys_tens[12 + site_idx, 0, 0] != 0.0:
          for att_idx, att in enumerate(self.attribute_list):
            self.crys_tens[att_idx, 12 + site_idx, :] *= self.normalize_near_max[att]
            self.crys_tens[12 + site_idx, att_idx, :] *= self.normalize_near_max[att]      

      for adj_idx in range(len(self.coord_list)):
        if site_idx != adj_idx:
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] *= self.normalize_near_max["length"]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] *= self.normalize_near_max["length"]

          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] *= self.normalize_near_max["dir"]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] -= self.direction_fix
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] *= self.normalize_near_max["dir"]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] -= self.direction_fix

    elif method == "avg_std":
      for site_idx in range(self.crys_tens.shape[0] - 12):

        #Check if an atom is present
        if self.crys_tens[12 + site_idx, 0, 0] != 0.0:
          for att_idx, att in enumerate(self.attribute_list):
            self.crys_tens[att_idx, 12 + site_idx, :] *= self.normalize_avg_std[att][2] + self.normalize_avg_std[att][3]
            self.crys_tens[att_idx, 12 + site_idx, :] -= self.normalize_avg_std[att][2]
            self.crys_tens[att_idx, 12 + site_idx, :] *= self.normalize_avg_std[att][1]
            self.crys_tens[att_idx, 12 + site_idx, :] += self.normalize_avg_std[att][0]
            
            self.crys_tens[12 + site_idx, att_idx, :] *= self.normalize_avg_std[att][2] + self.normalize_avg_std[att][3]
            self.crys_tens[12 + site_idx, att_idx, :] -= self.normalize_avg_std[att][2]
            self.crys_tens[12 + site_idx, att_idx, :] *= self.normalize_avg_std[att][1]
            self.crys_tens[12 + site_idx, att_idx, :] += self.normalize_avg_std[att][0] 

      for adj_idx in range(len(self.coord_list)):
        if site_idx != adj_idx:
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] *= self.normalize_avg_std["length"][2] + self.normalize_avg_std["length"][3]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] -= self.normalize_avg_std["length"][2]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] *= self.normalize_avg_std["length"][1]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 0] += self.normalize_avg_std["length"][0]
          
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] *= self.normalize_avg_std["length"][2] + self.normalize_avg_std["length"][3]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] -= self.normalize_avg_std["length"][2]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] *= self.normalize_avg_std["length"][1]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 0] += self.normalize_avg_std["length"][0]

          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] *= self.normalize_avg_std["dir"][2] + self.normalize_avg_std["dir"][3]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] -= self.normalize_avg_std["dir"][2]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] *= self.normalize_avg_std["dir"][1]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] += self.normalize_avg_std["dir"][0]
          self.crys_tens[12 + site_idx, 12 + adj_idx, 1:] -= self.direction_fix
          
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] *= self.normalize_avg_std["dir"][2] + self.normalize_avg_std["dir"][3]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] -= self.normalize_avg_std["dir"][2]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] *= self.normalize_avg_std["dir"][1]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] += self.normalize_avg_std["dir"][0]
          self.crys_tens[12 + adj_idx, 12 + site_idx, 1:] -= self.direction_fix

  def get_cif(self, file_path, ref_percent: float = 0.2, rel_coord_percent: float = None, num_coord_clusters: int = None, num_atom_clusters: int = 3, generator: str = "Real CIF"):
    crystal_cif = {}
    ref_angle = self.crys_tens[12, 7, 0]
    ref_avg = ref_angle
    for num in range(13, self.crys_tens.shape[0]):
      if abs(self.crys_tens[num, 7, 0] - ref_avg)/ref_avg >= ref_percent:
        break
      else:
        ref_angle += self.crys_tens[num, 7, 0]
        reference_average = ref_angle / (num - 11)
    
    constrained_crys_tens = self.crys_tens[:num, :num, :]
    for i in range(12, constrained_crys_tens.shape[0]):
      for j in range(12):
        avg_val = (np.sum(constrained_crys_tens[i, j, :]) + np.sum(constrained_crys_tens[j, i, :]))/(2 * constrained_crys_tens.shape[2]) 
        constrained_crys_tens[i, j, :], constrained_crys_tens[j, i, :] = avg_val, avg_val

    for i in range(4, 11):
      sum_val = 0
      for j in range(12, constrained_crys_tens.shape[0]):
        sum_val += np.sum(constrained_crys_tens[i, j, :]) + np.sum(constrained_crys_tens[j, i, :])
      avg_val = sum_val / ((constrained_crys_tens.shape[0] - 12) * constrained_crys_tens.shape[2])
      constrained_crys_tens[i, 12:, :], constrained_crys_tens[12:, i, :] = avg_val, avg_val
      crystal_cif[self.parameter_list[i - 4]] = avg_val

    for i in range(12, constrained_crys_tens.shape[0]):
      for j in range(12, constrained_crys_tens.shape[0]):
        avg_val = (constrained_crys_tens[i, j, 0] + constrained_crys_tens[j, i, 0])/2
        constrained_crys_tens[i, j, 0], constrained_crys_tens[j, i, 0] = avg_val, avg_val
        for k in range(1, 4):
          if constrained_crys_tens[i, j, k] > 0:
            avg_val = (constrained_crys_tens[i, j, k] + abs(constrained_crys_tens[j, i, k]))/2
            constrained_crys_tens[i, j, k], constrained_crys_tens[j, i, k] = avg_val, -avg_val
          else:
            avg_val = (abs(constrained_crys_tens[i, j, k]) + constrained_crys_tens[j, i, k])/2
            constrained_crys_tens[i, j, k], constrained_crys_tens[j, i, k] = -avg_val, avg_val

    crystal_cif["site_list"] = {}
    for i in range(12, constrained_crys_tens.shape[0]):
      crystal_cif["site_list"][i - 12] = {}
      crystal_cif["site_list"][i - 12]["atom"] = constrained_crys_tens[i, 0, 0]
      crystal_cif["site_list"][i - 12]["x"] = constrained_crys_tens[i, 1, 0]
      crystal_cif["site_list"][i - 12]["y"] = constrained_crys_tens[i, 2, 0]
      crystal_cif["site_list"][i - 12]["z"] = constrained_crys_tens[i, 3, 0]

      crystal_cif["site_list"][i - 12]["adj_list"] = []
      for j in range(12, constrained_crys_tens.shape[0]):
        adj_x = crystal_cif["site_list"][i - 12]["x"] - constrained_crys_tens[i, j, 1]
        adj_y = crystal_cif["site_list"][i - 12]["y"] - constrained_crys_tens[i, j, 2]
        adj_z = crystal_cif["site_list"][i - 12]["z"] - constrained_crys_tens[i, j, 3]
        crystal_cif["site_list"][i - 12]["adj_list"].append((adj_x, adj_y, adj_z))

    site_list = []
    atom_list = []
    if rel_coord_percent is None:
      rel_coord_percent = 1 - (1/constrained_crys_tens.shape[0])
    for site_idx in crystal_cif["site_list"]:
      site = crystal_cif["site_list"][site_idx]
      site_x = site["x"]
      site_y = site["y"]
      site_z = site["z"]
      adj_x_list = []
      adj_y_list = []
      adj_z_list = []
      adj_coord_list =site["adj_list"]
      for adj_idx in range(len(adj_coord_list)):
        adj_coord = adj_coord_list[adj_idx]
        adj_x_list.append(adj_coord[0])
        adj_y_list.append(adj_coord[1])
        adj_z_list.append(adj_coord[2])
      
      site_coord_percent = 1 - rel_coord_percent
      x_rel = site_coord_percent*site_x + rel_coord_percent*np.average(adj_x_list)
      y_rel = site_coord_percent*site_y + rel_coord_percent*np.average(adj_y_list)
      z_rel = site_coord_percent*site_z + rel_coord_percent*np.average(adj_z_list)
      atom_list.append(np.around(site["atom"]))
      site_list.append((x_rel, y_rel, z_rel))

    kmeans_coord_list = []
    for i in range(len(site_list)):
      for j in range(3):
        kmeans_coord_list.append(site_list[i][j])

    if num_coord_clusters is None:
      num_coord_clusters = len(kmeans_coord_list)

    kmeans_coord = KMeans(n_clusters = num_coord_clusters).fit(np.array(kmeans_coord_list).reshape(-1, 1))
    for i in range(len(kmeans_coord_list)):
      kmeans_coord_list[i] = kmeans_coord.cluster_centers_[kmeans_coord.labels_[i]][0]
    
    for i in range(0, len(kmeans_coord_list), 3):
      site_list[i//3] = [kmeans_coord_list[i], kmeans_coord_list[i + 1], kmeans_coord_list[i + 2]]

    #TODO Finish PotScoring

    num_atom_clusters = min(num_atom_clusters, len(atom_list))
    kmeans_atom = KMeans(n_clusters=int(num_atom_clusters)).fit(np.array(atom_list).reshape(-1, 1))

    for i in range(len(atom_list)):
      atom_list[i] = np.around(kmeans_atom.cluster_centers_[kmeans_atom.labels_[i]])[0]

    lattice = Lattice.from_parameters(a = crystal_cif["a"], b = crystal_cif["b"], c = crystal_cif["c"], alpha = crystal_cif["alpha"], beta = crystal_cif["beta"], gamma = crystal_cif["gamma"])
    struct = Structure(lattice = lattice, species = atom_list, coords = site_list, to_unit_cell=True)
    written_cif = str(CifWriter(struct))
    with open(file_path, "w") as file:
      file.write("Generated by: " + generator + "\n" + "Num unique sites: " + str(num_coord_clusters) + "\n" + "Num unique elements: " + str(num_atom_clusters) + "\n\n" + written_cif)

  def construct_crys_tens(self):
    self.crys_tens = np.zeros((64, 64, 4))
    for site_idx in range(len(self.coord_list)):
      #Insert atomic number of atom
      self.crys_tens[0, 12 + site_idx, :] = self.element_list[site_idx]
      self.crys_tens[12 + site_idx, 0, :] = self.element_list[site_idx]

      #Insert fractional coordinates of atom (x, y, z)
      self.crys_tens[1, 12 + site_idx, :] = self.coord_list[site_idx][0]
      self.crys_tens[12 + site_idx, 1, :] = self.coord_list[site_idx][0]
      self.crys_tens[2, 12 + site_idx, :] = self.coord_list[site_idx][1]
      self.crys_tens[12 + site_idx, 2, :] = self.coord_list[site_idx][1]
      self.crys_tens[3, 12 + site_idx, :] = self.coord_list[site_idx][2]
      self.crys_tens[12 + site_idx, 3, :] = self.coord_list[site_idx][2]

      #Insert lattice parameters (a, b, c)
      self.crys_tens[4, 12 + site_idx, :] = self.a
      self.crys_tens[12 + site_idx, 4, :] = self.a
      self.crys_tens[5, 12 + site_idx, :] = self.b
      self.crys_tens[12 + site_idx, 5, :] = self.b
      self.crys_tens[6, 12 + site_idx, :] = self.c
      self.crys_tens[12 + site_idx, 6, :] = self.c

      #Insert lattice parameters (alpha, beta, gamma)
      self.crys_tens[7, 12 + site_idx, :] = self.alpha
      self.crys_tens[12 + site_idx, 7, :] = self.alpha
      self.crys_tens[8, 12 + site_idx, :] = self.beta
      self.crys_tens[12 + site_idx, 8, :] = self.beta
      self.crys_tens[9, 12 + site_idx, :] = self.gamma
      self.crys_tens[12 + site_idx, 9, :] = self.gamma

      #Insert space group number
      self.crys_tens[10, 12 + site_idx, :] = self.sg
      self.crys_tens[12 + site_idx, 10, :] = self.sg

      for adj_idx in range(len(self.coord_list)):
        x_diff = self.coord_list[site_idx][0] - self.coord_list[adj_idx][0]
        y_diff = self.coord_list[site_idx][1] - self.coord_list[adj_idx][1]
        z_diff = self.coord_list[site_idx][2] - self.coord_list[adj_idx][2]
        length = math.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

        #Insert adjacency matrix on first layer
        self.crys_tens[12 + site_idx, 12 + adj_idx, 0] = length
        self.crys_tens[12 + adj_idx, 12 + site_idx, 0] = length

        #Insert the dimensional differences in the latter three layers (x, y, z)
        self.crys_tens[12 + site_idx, 12 + adj_idx, 1] = x_diff
        self.crys_tens[12 + adj_idx, 12 + site_idx, 1] = -x_diff
        self.crys_tens[12 + site_idx, 12 + adj_idx, 2] = y_diff
        self.crys_tens[12 + adj_idx, 12 + site_idx, 2] = -y_diff
        self.crys_tens[12 + site_idx, 12 + adj_idx, 3] = z_diff
        self.crys_tens[12 + adj_idx, 12 + site_idx, 3] = -z_diff

  def get_crys_tens(self, normalized, method):
    if normalized:
      self.normalize_crys_tens(method)
      return self.crys_tens
    else:
      self.unnormalize_crys_tens(method)
      return self.crys_tens


class StackedCrysTensor:
  def __init__(self, *args):
    if len(args) == 0:
      raise ValueError("Please input a directory of CIF files or a list of CrysTens'.")
    self.crys_tensor = []
    if isinstance(args[0], str):
      for path in os.listdir(args[0]):
        cif_path = os.path.join(args[0], path)
        self.crys_tensor.append(CrysTens(cif_path))
    else:
      self.crys_tensor = args[0]
  
  def get_stacked_crys_tensor(self, normalized: str = True, method: str = "avg_std"):
    crys_tensor_np = np.zeros((len(self.crys_tensor), 64, 64, 4))
    for idx, crys_tens in enumerate(self.crys_tensor):
      crys_tensor_np[idx, :, :, :] = crys_tens.get_crys_tens(normalized, method)
    return crys_tensor_np