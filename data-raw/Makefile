# file directories
DATA=../data

# all datasets for package
all: $(DATA)/therapy-pars.rda $(DATA)/pain.rda $(DATA)/util-pars.rda $(DATA)/duration-eular-pars.rda
all: $(DATA)/acr2eular.rda $(DATA)/eular2haq.rda $(DATA)/acr2haq.rda $(DATA)/haq-lprog-age.rda 
all: $(DATA)/mort-or.rda $(DATA)/mort-hr-haqdif.rda $(DATA)/lifetable-male.rda $(DATA)/lifetable-female.rda 
all: $(DATA)/mgmt-cost.rda

# parameters that vary by treatment
TP_DEPS := therapy-pars.R func.R cost.xlsx 
TP_TARGET := $(DATA)/therapy-pars.rda figs/infections-grijalva-cumhaz.pdf figs/infections-van-dartel-cumhaz.pdf
$(TP_TARGET): $(TP_DEPS)
	Rscript therapy-pars.R

# utility
$(DATA)/pain.rda $(DATA)/util-pars.rda $(DATA)/util-wailoo-pars.rda: utility.R coef-alava-2013.csv vcov-alava-2013.csv
	Rscript utility.R 

# treatment duration
TD_DEPS := duration-eular.R func.R duration-bsrbr-eular-good-ggamma.csv duration-bsrbr-eular-good-risk.csv 
TD_DEPS += duration-bsrbr-eular-moderate-ggamma.csv duration-bsrbr-eular-moderate-risk.csv
TD_DEPS += duration-corrona-allpat.csv duration-corrona-allpat-risk.csv
TD_TARGET := $(DATA)/duration-eular-pars.rda 
$(TD_TARGET): $(TD_DEPS)
	Rscript duration-eular.R
		
# Mappings
$(DATA)/acr2haq.rda $(DATA)/acr2eular.rda: mappings.R
	Rscript mappings.R
	
# HAQ change by EULAR response
$(DATA)/eular2haq.rda: eular2haq.R
	Rscript eular2haq.R
	
# HAQ progression rates by age
$(DATA)/haq-lprog-age.rda: haq-lprog-age.R
	Rscript haq-lprog-age.R
	
# Mortality odds ratios
$(DATA)/mort-or.rda: mort-or.R
	Rscript mort-or.R
	
# Mortality hazard ratios for change in HAQ
$(DATA)/mort-hr-haqdif.rda: mort-hr-haqdif.R
	Rscript mort-hr-haqdif.R
	
# Lifetables
$(DATA)/lifetable-male.rda: lifetable.R
	Rscript lifetable.R
$(DATA)/lifetable-female.rda: lifetable.R
	Rscript lifetable.R
	
# Management costs
$(DATA)/mgmt-cost.rda: mgmt-cost.R
	Rscript mgmt-cost.R
	