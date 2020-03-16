#This is a Dockerfile used to create a containter for testing used by test action.
FROM continuumio/miniconda3


RUN conda install -c conda-forge -c intbio -c mdanalysis pandas pytest


CMD [ "/bin/bash" ]