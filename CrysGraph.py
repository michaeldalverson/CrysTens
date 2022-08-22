class CrysGraph:
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
    
    self.construct_crysgraph()

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

    self.construct_crysgraph()

  def from_CIF_obj(self):
    #TODO
    pass

  def normalize_crysgraph(self, method: str = "divide_near_max"):
    if self.normalized:
      return
    self.normalized = True
    if method == "divide_near_max":
      print(self.crysgraph.shape)
      for site_idx in range(self.crysgraph.shape[0] - 12):

        #Check if an atom is present
        if self.crysgraph[12 + site_idx, 0, 0] != 0.0:
          self.crysgraph[0, 12 + site_idx, :] /= self.normalize_near_max["atom"]
          self.crysgraph[12 + site_idx, 0, :] /= self.normalize_near_max["atom"]

          self.crysgraph[1, 12 + site_idx, :] /= self.normalize_near_max["x"]
          self.crysgraph[12 + site_idx, 1, :] /= self.normalize_near_max["x"]
          self.crysgraph[2, 12 + site_idx, :] /= self.normalize_near_max["y"]
          self.crysgraph[12 + site_idx, 2, :] /= self.normalize_near_max["y"]
          self.crysgraph[3, 12 + site_idx, :] /= self.normalize_near_max["z"]
          self.crysgraph[12 + site_idx, 3, :] /= self.normalize_near_max["z"]

          self.crysgraph[4, 12 + site_idx, :] /= self.normalize_near_max["a"]
          self.crysgraph[12 + site_idx, 4, :] /= self.normalize_near_max["a"]
          self.crysgraph[5, 12 + site_idx, :] /= self.normalize_near_max["b"]
          self.crysgraph[12 + site_idx, 5, :] /= self.normalize_near_max["b"]
          self.crysgraph[6, 12 + site_idx, :] /= self.normalize_near_max["c"]
          self.crysgraph[12 + site_idx, 6, :] /= self.normalize_near_max["c"]

          self.crysgraph[7, 12 + site_idx, :] /= self.normalize_near_max["alpha"]
          self.crysgraph[12 + site_idx, 7, :] /= self.normalize_near_max["alpha"]
          self.crysgraph[8, 12 + site_idx, :] /= self.normalize_near_max["beta"]
          self.crysgraph[12 + site_idx, 8, :] /= self.normalize_near_max["beta"]
          self.crysgraph[9, 12 + site_idx, :] /= self.normalize_near_max["gamma"]
          self.crysgraph[12 + site_idx, 9, :] /= self.normalize_near_max["gamma"]   

          self.crysgraph[12 + site_idx, 10, :] /= self.normalize_near_max["sg"]        

      for adj_idx in range(len(self.coord_list)):

        if site_idx != adj_idx:
          self.crysgraph[12 + site_idx, 12 + adj_idx, 0] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 0] /= self.normalize_near_max["dir"]

          self.crysgraph[12 + site_idx, 12 + adj_idx, 1] += self.direction_fix
          self.crysgraph[12 + site_idx, 12 + adj_idx, 1] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 1] += self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 1] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + site_idx, 12 + adj_idx, 2] += self.direction_fix
          self.crysgraph[12 + site_idx, 12 + adj_idx, 2] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 2] += self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 2] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + site_idx, 12 + adj_idx, 3] += self.direction_fix
          self.crysgraph[12 + site_idx, 12 + adj_idx, 3] /= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 3] += self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 3] /= self.normalize_near_max["dir"]

  def unnormalize_crysgraph(self, method: str = "divide_near_max"):
    if not self.normalized:
      return
    self.normalized = False
    if method == "divide_near_max":
      for site_idx in range(self.crysgraph.shape[0]):

        #Check if an atom is present
        print(self.crysgraph)
        if self.crysgraph[12 + site_idx, 0, 0] != 0.0:
          self.crysgraph[0, 12 + site_idx, :] *= self.normalize_near_max["atom"]
          self.crysgraph[12 + site_idx, 0, :] *= self.normalize_near_max["atom"]

          self.crysgraph[1, 12 + site_idx, :] *= self.normalize_near_max["x"]
          self.crysgraph[12 + site_idx, 1, :] *= self.normalize_near_max["x"]
          self.crysgraph[2, 12 + site_idx, :] *= self.normalize_near_max["y"]
          self.crysgraph[12 + site_idx, 2, :] *= self.normalize_near_max["y"]
          self.crysgraph[3, 12 + site_idx, :] *= self.normalize_near_max["z"]
          self.crysgraph[12 + site_idx, 3, :] *= self.normalize_near_max["z"]

          self.crysgraph[4, 12 + site_idx, :] *= self.normalize_near_max["a"]
          self.crysgraph[12 + site_idx, 4, :] *= self.normalize_near_max["a"]
          self.crysgraph[5, 12 + site_idx, :] *= self.normalize_near_max["b"]
          self.crysgraph[12 + site_idx, 5, :] *= self.normalize_near_max["b"]
          self.crysgraph[6, 12 + site_idx, :] *= self.normalize_near_max["c"]
          self.crysgraph[12 + site_idx, 6, :] *= self.normalize_near_max["c"]

          self.crysgraph[7, 12 + site_idx, :] *= self.normalize_near_max["alpha"]
          self.crysgraph[12 + site_idx, 7, :] *= self.normalize_near_max["alpha"]
          self.crysgraph[8, 12 + site_idx, :] *= self.normalize_near_max["beta"]
          self.crysgraph[12 + site_idx, 8, :] *= self.normalize_near_max["beta"]
          self.crysgraph[9, 12 + site_idx, :] *= self.normalize_near_max["gamma"]
          self.crysgraph[12 + site_idx, 9, :] *= self.normalize_near_max["gamma"]   

          self.crysgraph[12 + site_idx, 10, :] *= self.normalize_near_max["sg"]        

      for adj_idx in range(len(self.coord_list)):

        if site_idx != adj_idx:
          self.crysgraph[12 + site_idx, 12 + adj_idx, 0] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 0] *= self.normalize_near_max["dir"]

          self.crysgraph[12 + site_idx, 12 + adj_idx, 1] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + site_idx, 12 + adj_idx, 1] -= self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 1] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 1] -= self.direction_fix
          self.crysgraph[12 + site_idx, 12 + adj_idx, 2] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + site_idx, 12 + adj_idx, 2] -= self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 2] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 2] -= self.direction_fix
          self.crysgraph[12 + site_idx, 12 + adj_idx, 3] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + site_idx, 12 + adj_idx, 3] -= self.direction_fix
          self.crysgraph[12 + adj_idx, 12 + site_idx, 3] *= self.normalize_near_max["dir"]
          self.crysgraph[12 + adj_idx, 12 + site_idx, 3] -= self.direction_fix

  def construct_crysgraph(self):
    self.crysgraph = np.zeros((64, 64, 4))
    for site_idx in range(len(self.coord_list)):
      #Insert atomic number of atom
      self.crysgraph[0, 12 + site_idx, :] = self.element_list[site_idx]
      self.crysgraph[12 + site_idx, 0, :] = self.element_list[site_idx]

      #Insert fractional coordinates of atom (x, y, z)
      self.crysgraph[1, 12 + site_idx, :] = self.coord_list[site_idx][0]
      self.crysgraph[12 + site_idx, 1, :] = self.coord_list[site_idx][0]
      self.crysgraph[2, 12 + site_idx, :] = self.coord_list[site_idx][1]
      self.crysgraph[12 + site_idx, 2, :] = self.coord_list[site_idx][1]
      self.crysgraph[3, 12 + site_idx, :] = self.coord_list[site_idx][2]
      self.crysgraph[12 + site_idx, 3, :] = self.coord_list[site_idx][2]

      #Insert lattice parameters (a, b, c)
      self.crysgraph[4, 12 + site_idx, :] = self.a
      self.crysgraph[12 + site_idx, 4, :] = self.a
      self.crysgraph[5, 12 + site_idx, :] = self.b
      self.crysgraph[12 + site_idx, 5, :] = self.b
      self.crysgraph[6, 12 + site_idx, :] = self.c
      self.crysgraph[12 + site_idx, 6, :] = self.c

      #Insert lattice parameters (alpha, beta, gamma)
      self.crysgraph[7, 12 + site_idx, :] = self.alpha
      self.crysgraph[12 + site_idx, 7, :] = self.alpha
      self.crysgraph[8, 12 + site_idx, :] = self.beta
      self.crysgraph[12 + site_idx, 8, :] = self.beta
      self.crysgraph[9, 12 + site_idx, :] = self.gamma
      self.crysgraph[12 + site_idx, 9, :] = self.gamma

      #Insert space group number
      self.crysgraph[10, 12 + site_idx, :] = self.sg
      self.crysgraph[12 + site_idx, 10, :] = self.sg

      for adj_idx in range(len(self.coord_list)):
        x_diff = self.coord_list[site_idx][0] - self.coord_list[adj_idx][0]
        y_diff = self.coord_list[site_idx][1] - self.coord_list[adj_idx][1]
        z_diff = self.coord_list[site_idx][2] - self.coord_list[adj_idx][2]
        length = math.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

        #Insert adjacency matrix on first layer
        self.crysgraph[12 + site_idx, 12 + adj_idx, 0] = length
        self.crysgraph[12 + adj_idx, 12 + site_idx, 0] = length

        #Insert the dimensional differences in the latter three layers (x, y, z)
        self.crysgraph[12 + site_idx, 12 + adj_idx, 1] = x_diff
        self.crysgraph[12 + adj_idx, 12 + site_idx, 1] = x_diff
        self.crysgraph[12 + site_idx, 12 + adj_idx, 2] = y_diff
        self.crysgraph[12 + adj_idx, 12 + site_idx, 2] = y_diff
        self.crysgraph[12 + site_idx, 12 + adj_idx, 3] = z_diff
        self.crysgraph[12 + adj_idx, 12 + site_idx, 3] = z_diff

  def get_crys_graph(self, normalized):
    if normalized:
      self.normalize_crysgraph()
      return self.crysgraph
    else:
      self.unnormalize_crysgraph()
      return self.crysgraph


class CrysTensor:
  def __init__(self, *args):
    if len(args) == 0:
      raise ValueError("Please input a directory of CIF files or a list of CrysGraphs.")
    self.crys_tensor = []
    if isinstance(args[0], str):
      for path in os.listdir(args[0]):
        cif_path = os.path.join(args[0], path)
        print(cif_path)
        self.crys_tensor.append(CrysGraph(cif_path))
    else:
      self.crys_tensor = args[0]
  
  def get_crys_tensor(self, normalized: str = True):
    crys_tensor_np = np.zeros((len(self.crys_tensor), 64, 64, 4))
    for idx, crys_graph in enumerate(self.crys_tensor):
      print(crys_graph)
      crys_tensor_np[idx, :, :, :] = crys_graph.get_crys_graph(normalized)
    return crys_tensor_np