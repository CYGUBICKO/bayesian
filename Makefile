## Target
current: target
-include target.mk

## Makestuff setup
Sources += Makefile 
msrepo = https://github.com/dushoff
ms = makestuff
-include $(ms)/os.mk

Ignore += $(ms)
Makefile: $(ms)/Makefile
$(ms)/Makefile:
	git clone $(msrepo)/$(ms)
	ls $@

######################################################################

Sources += $(wildcard *.R *.md)

######################################################################

## These are codes based on the Bayesian Nonparametric data analysis

# Example 3: T-Cell Receptors Guindani et al.

tcell_ex3.Rout: tcell_ex3.R

# Example 4: Gene expression data
## ln -fs ~/grive/McMaster/PhD/Readings-Materials/Books/Book-Codes/Nonparametric-bayesian
## ln -fs ~/grive/McMaster/PhD/Readings-Materials/Books/Book-Codes/Nonparametric-bayesian/Epildata.txt

gene_ex4.Rout: gene_ex4.R

epilepsy.Rout: epilepsy.R

## MCMC
simpleMCMC.Rout: simpleMCMC.R

myMCMC.Rout: myMCMC.R


######################################################################


### Makestuff rules
-include $(ms)/git.mk
-include $(ms)/git.mk
-include $(ms)/visual.mk
-include $(ms)/wrapR.mk
-include $(ms)/texdeps.mk
-include $(ms)/pandoc.mk
-include $(ms)/autorefs.mk

