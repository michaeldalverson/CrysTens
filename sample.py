from CrysGraph import CrysTensor as CG

graph_generator = CG

t = graph_generator("test\sample_cifs").get_crys_tensor

print(t)
