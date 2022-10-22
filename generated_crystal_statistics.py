import numpy as np
from typing import List
from sklearn.cluster import KMeans
from pymatgen.io.cif import CifWriter
from pymatgen.core import Lattice, Structure

direction_fix = 0.9999933
normalize_avg_std = {'atom': [22.417031299395948, 27.503517675002012, 0.778701530199605, 2.7844790475732646],
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

attribute_list = ["atom", "x", "y", "z", "a", "b", "c", "alpha", "beta", "gamma", "sg"]
parameter_list = ["a", "b", "c", "alpha", "beta", "gamma", "sg"]

def unnormalize_crysgraph(crysgraph):
  for site_idx in range(crysgraph.shape[0] - 12):
    #Check if an atom is present
    if crysgraph[12 + site_idx, 0, 0] != 0.0:
      for att_idx, att in enumerate(attribute_list):
        crysgraph[att_idx, 12 + site_idx, :] *= normalize_avg_std[att][2] + normalize_avg_std[att][3]
        crysgraph[att_idx, 12 + site_idx, :] -= normalize_avg_std[att][2]
        crysgraph[att_idx, 12 + site_idx, :] *= normalize_avg_std[att][1]
        crysgraph[att_idx, 12 + site_idx, :] += normalize_avg_std[att][0]
        
        crysgraph[12 + site_idx, att_idx, :] *= normalize_avg_std[att][2] + normalize_avg_std[att][3]
        crysgraph[12 + site_idx, att_idx, :] -= normalize_avg_std[att][2]
        crysgraph[12 + site_idx, att_idx, :] *= normalize_avg_std[att][1]
        crysgraph[12 + site_idx, att_idx, :] += normalize_avg_std[att][0] 

  for adj_idx in range(crysgraph.shape[0] - 12):
    if site_idx != adj_idx:
      crysgraph[12 + site_idx, 12 + adj_idx, 0] *= normalize_avg_std["length"][2] + normalize_avg_std["length"][3]
      crysgraph[12 + site_idx, 12 + adj_idx, 0] -= normalize_avg_std["length"][2]
      crysgraph[12 + site_idx, 12 + adj_idx, 0] *= normalize_avg_std["length"][1]
      crysgraph[12 + site_idx, 12 + adj_idx, 0] += normalize_avg_std["length"][0]
      
      crysgraph[12 + adj_idx, 12 + site_idx, 0] *= normalize_avg_std["length"][2] + normalize_avg_std["length"][3]
      crysgraph[12 + adj_idx, 12 + site_idx, 0] -= normalize_avg_std["length"][2]
      crysgraph[12 + adj_idx, 12 + site_idx, 0] *= normalize_avg_std["length"][1]
      crysgraph[12 + adj_idx, 12 + site_idx, 0] += normalize_avg_std["length"][0]

      crysgraph[12 + site_idx, 12 + adj_idx, 1:] *= normalize_avg_std["dir"][2] + normalize_avg_std["dir"][3]
      crysgraph[12 + site_idx, 12 + adj_idx, 1:] -= normalize_avg_std["dir"][2]
      crysgraph[12 + site_idx, 12 + adj_idx, 1:] *= normalize_avg_std["dir"][1]
      crysgraph[12 + site_idx, 12 + adj_idx, 1:] += normalize_avg_std["dir"][0]
      crysgraph[12 + site_idx, 12 + adj_idx, 1:] -= direction_fix
      
      crysgraph[12 + adj_idx, 12 + site_idx, 1:] *= normalize_avg_std["dir"][2] + normalize_avg_std["dir"][3]
      crysgraph[12 + adj_idx, 12 + site_idx, 1:] -= normalize_avg_std["dir"][2]
      crysgraph[12 + adj_idx, 12 + site_idx, 1:] *= normalize_avg_std["dir"][1]
      crysgraph[12 + adj_idx, 12 + site_idx, 1:] += normalize_avg_std["dir"][0]
      crysgraph[12 + adj_idx, 12 + site_idx, 1:] -= direction_fix
  return crysgraph

  def get_cif(self, file_path, ref_percent: float = 0.2, rel_coord_percent: float = None, num_coord_clusters: int = None, num_atom_clusters: int = 3, generator: str = "Real CIF", dir_diff_avg: List = None):
    crystal_cif = {}
    ref_angle = self.crys_graph[12, 7, 0]
    ref_avg = ref_angle
    for num in range(13, self.crys_graph.shape[0]):
      if abs(self.crys_graph[num, 7, 0] - ref_avg)/ref_avg >= ref_percent:
        break
      else:
        ref_angle += self.crys_graph[num, 7, 0]
        reference_average = ref_angle / (num - 11)
    
    constrained_crysgraph = self.crysgraph[:num, :num, :]
    for i in range(12, constrained_crysgraph.shape[0]):
      for j in range(12):
        avg_val = (np.sum(constrained_crysgraph[i, j, :]) + np.sum(constrained_crysgraph[j, i, :]))/(2 * constrained_crysgraph.shape[2]) 
        constrained_crysgraph[i, j, :], constrained_crysgraph[j, i, :] = avg_val, avg_val

    for i in range(4, 11):
      sum_val = 0
      for j in range(12, constrained_crysgraph.shape[0]):
        sum_val += np.sum(constrained_crysgraph[i, j, :]) + np.sum(constrained_crysgraph[j, i, :])
      avg_val = sum_val / ((constrained_crysgraph.shape[0] - 12) * constrained_crysgraph.shape[2])
      constrained_crysgraph[i, 12:, :], constrained_crysgraph[12:, i, :] = avg_val, avg_val
      crystal_cif[self.parameter_list[i - 4]] = avg_val

    for i in range(12, constrained_crysgraph.shape[0]):
      for j in range(12, constrained_crysgraph.shape[0]):
        avg_val = (constrained_crysgraph[i, j, 0] + constrained_crysgraph[j, i, 0])/2
        constrained_crysgraph[i, j, 0], constrained_crysgraph[j, i, 0] = avg_val, avg_val
        for k in range(1, 4):
          if constrained_crysgraph[i, j, k] > 0:
            avg_val = (constrained_crysgraph[i, j, k] + abs(constrained_crysgraph[j, i, k]))/2
            constrained_crysgraph[i, j, k], constrained_crysgraph[j, i, k] = avg_val, -avg_val
          else:
            avg_val = (abs(constrained_crysgraph[i, j, k]) + constrained_crysgraph[j, i, k])/2
            constrained_crysgraph[i, j, k], constrained_crysgraph[j, i, k] = -avg_val, avg_val

    crystal_cif["site_list"] = {}
    for i in range(12, constrained_crysgraph.shape[0]):
      crystal_cif["site_list"][i - 12] = {}
      crystal_cif["site_list"][i - 12]["atom"] = constrained_crysgraph[i, 0, 0]
      crystal_cif["site_list"][i - 12]["x"] = constrained_crysgraph[i, 1, 0]
      crystal_cif["site_list"][i - 12]["y"] = constrained_crysgraph[i, 2, 0]
      crystal_cif["site_list"][i - 12]["z"] = constrained_crysgraph[i, 3, 0]

      crystal_cif["site_list"][i - 12]["adj_list"] = []
      for j in range(12, constrained_crysgraph.shape[0]):
        adj_x = crystal_cif["site_list"][i - 12]["x"] - constrained_crysgraph[i, j, 1]
        adj_y = crystal_cif["site_list"][i - 12]["y"] - constrained_crysgraph[i, j, 2]
        adj_z = crystal_cif["site_list"][i - 12]["z"] - constrained_crysgraph[i, j, 3]
        crystal_cif["site_list"][i - 12]["adj_list"].append((adj_x, adj_y, adj_z))

    site_list = []
    atom_list = []
    if rel_coord_percent is None:
      rel_coord_percent = 1 - (1/constrained_crysgraph.shape[0])
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
      if dir_diff_avg is not None:
        dir_diff_avg.append(abs(site_x - x_rel)/adj_x_list)
        dir_diff_avg.append(abs(site_y - y_rel)/adj_y_list)
        dir_diff_avg.append(abs(site_z - z_rel)/adj_z_list)

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

def get_statistics(crysgraph_path, cif_folder):
    xyzvar = []
    paramvar = []
    anglevar = []
    sgvar = []
    dir_diff_avg = []
    for crys in (np_crys):
    if np_crys.shape != (64, 64, 4):
        new_crys = np.zeros((64, 64, 4))
        for i in range(4):
        new_crys[:, :, i] = crys[i, :, :]
        crys = new_crys[:]
    crysgraph = unnormalize_crysgraph(crys)

    get_cif(file_path, 0.2, None, None, 3, "Real CIF", dir_diff_avg)
    try:
        load_crystal = Structure.from_file(file_path)
    except:
        continue

    ref_angle = crysgraph[12, 7, 0]
    ref_avg = ref_angle
    for num in range(13, crysgraph.shape[0]):
        if abs(crysgraph[num, 7, 0] - ref_avg)/ref_avg >= 0.2:
        break
        else:
        ref_angle += crysgraph[num, 7, 0]
        reference_average = ref_angle / (num - 11)
    crysgraph = crysgraph[:num, :num, :]

    for i in range(12, crysgraph.shape[0]):
        for j in range(1, 3):
        xyz = []
        for k in range(4):
            xyz.append(crysgraph[i, j, k])
            xyz.append(crysgraph[j, i, k])
        xyzvar.append(np.var(xyz))
    
    for i in range(4, 11):
        row = crysgraph[i, 12:, :]
        col = crysgraph[12:, i, :]
        total = []
        for j in range(4):
        total += row[j, :]
        total += col[j, :]
        if 3 < i < 7:
            paramvar.append(np.var(total))
        elif 7 <= i < 10:
            anglevar.append(np.var(total))
        else:
            sgvar.append(np.var(total))

