# CrysTens
CrysTens is a representation for storing crystal structure information that is originally in the form of Crystallographic Information Files (CIFs). CrysTens is a tensor of size 64x64x4 that can be used in any type of machine learning application involving crystal structures. This repository houses code for creating a stack of CrysTens', using the stack to train either a Vanilla Generative Adversarial Network (GAN), a Wasserstein GAN, or a diffusion model from https://github.com/lucidrains/imagen-pytorch/tree/main/imagen_pytorch. Once a model has been trained, newly generated CIFs can be created from a stack of generated CrysTens'. The details of the CrysTens representation, model training, and model analysis in the field of material discovery can be found here: https://chemrxiv.org/engage/chemrxiv/article-details/63694d64fbfd387c25d2d395.

## Getting Started
### Making a Stacked CrysTens
In order to train a CrysTens generative model for material discovery, a stack of CrysTens' is required. Using ```get_stacked_crys_tens.py``` and a Crystal Dictionary, any size of Stacked CrysTens can be received. 

#### Crystal Dictionary