# file directories
DATA=../data

# all datasets for package
all: $(DATA)/therapy-pars.rda $(DATA)/pain.rda $(DATA)/util-mixture-pars.rda $(DATA)/ttd-pars.rda
all: $(DATA)/acr2eular.rda $(DATA)/eular2haq.rda $(DATA)/acr2haq.rda $(DATA)/haq-lprog-age.rda 
all: $(DATA)/mort-or.rda $(DATA)/mort-hr-haqdif.rda $(DATA)/lifetable-male.rda $(DATA)/lifetable-female.rda 
all: $(DATA)/mgmt-cost.rda $(DATA)/haq-lcgm-pars.rda

# parameters that vary by treatment
TP_DEPS := therapy-pars.R func.R cost.xlsx 
TP_TARGET := $(DATA)/therapy-pars.rda figs/infections-grijalva-cumhaz.pdf figs/infections-van-dartel-cumhaz.pdf
$(TP_TARGET): $(TP_DEPS)
	Rscript therapy-pars.R

# HAQ LCGM
$(DATA)/haq-lcgm-pars.rda: haq-lcgm.R coef-norton-2014.csv
	Rscript haq-lcgm.R

# utility
$(DATA)/pain.rda $(DATA)/util-mixture-pars.rda $(DATA)/util-wailoo-pars.rda: utility.R coef-alava-2013.csv vcov-alava-2013.csv
	Rscript utility.R 

# treatment duration
TD_DEPS := ttd.R func.R ttd-bsrbr-eular-good-ggamma.csv ttd-bsrbr-eular-good-risk.csv 
TD_DEPS += ttd-bsrbr-eular-moderate-ggamma.csv ttd-bsrbr-eular-moderate-risk.csv
TD_DEPS += ttd-corrona-allpat.csv ttd-corrona-allpat-risk.csv
TD_TARGET := $(DATA)/ttd-pars.rda 
$(TD_TARGET): $(TD_DEPS)
	Rscript ttd.R
		
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
	