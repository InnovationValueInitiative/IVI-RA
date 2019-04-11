# file directories
DATA=../data

# All datasets, tables, and figures for package
all: $(DATA)/treatments.rda 
all: $(DATA)/nma-acr-naive.rda $(DATA)/nma-acr-naive-nice.rda 
all: $(DATA)/nma-das28-naive.rda $(DATA)/nma-haq-naive.rda
all: $(DATA)/acr2das28.rda $(DATA)/acr2sdai.rda $(DATA)/acr2cdai.rda 
all: $(DATA)/acr2eular.rda $(DATA)/eular2haq.rda $(DATA)/acr2haq.rda 
all: $(DATA)/haq-lprog.rda  $(DATA)/haq-lcgm.rda
all: $(DATA)/ttd.rda $(DATA)/ttsi.rda
all: $(DATA)/mort-or.rda $(DATA)/mort-hr-haqdif.rda $(DATA)/lifetable-male.rda $(DATA)/lifetable-female.rda 
all: $(DATA)/tx-cost.rda $(DATA)/hosp-cost.rda $(DATA)/mgmt-cost.rda $(DATA)/prod-loss.rda 
all: $(DATA)/utility-tx-attr.rda $(DATA)/utility-mixture.rda $(DATA)/utility-wailoo.rda
all: $(DATA)/incidence-female.rda $(DATA)/incidence-male.rda 
all: figs/ttd-corrona-ipdsurv-gengamma.pdf 
all: figs/ttd-bsrbr-ipdsurv-comp-gengamma-by-eular.pdf
all: figs/haq-lcgm-obsexp.pdf
all: tables/ttd-corrona-ipdsurv-ic.csv tables/ttd-bsrbr-ipdsurv-eular-ic.csv
all: tables/ttd-bsrbr-ipdsurv-eular-adjusted-ic.csv

# Treatment lookup
$(DATA)/treatments.rda: treatments.R treatments.csv
	Rscript treatments.R
	
# NMA
NMA_DEPS := treatments.csv nma.R nma-acr-naive-re-coda.csv nma-acr-naive-re-crosswalk.csv
NMA_DEPS += nma-das28-naive-re-coda.csv nma-das28-naive-re-crosswalk.csv
NMA_DEPS += nma-haq-naive-re-coda.csv nma-haq-naive-re-crosswalk.csv
NMA_TARGET := $(DATA)/nma-acr-naive.rda $(DATA)/nma-acr-naive-nice.rda 
NMA_TARGET +=  $(DATA)/nma-das28-naive.rda $(DATA)/nma-haq-naive.rda
$(NMA_TARGET): $(NMA_DEPS)
	Rscript nma.R

# Mappings
MAP_TARGET := := $(DATA)/acr2das28.rda $(DATA)/acr2sdai.rda $(DATA)/acr2cdai.rda 
MAP_TARGET += $(DATA)/acr2eular.rda $(DATA)/eular2haq.rda $(DATA)/acr2haq.rda
$(MAP_TARGET): mappings.R
	Rscript mappings.R

# HAQ progression rates
$(DATA)/haq-lprog.rda: haq-lprog.R treatments.csv
	Rscript haq-lprog.R
	
# HAQ LCGM
$(DATA)/haq-lcgm.rda figs/haq-lcgm-obsexp.pdf: haq-lcgm.R coef-norton-2014.csv
	Rscript haq-lcgm.R
	
# Treatment duration
TD_DEPS := ttd.R func.R ttd-bsrbr-eular-good-ggamma.csv ttd-bsrbr-eular-good-risk.csv 
TD_DEPS += ttd-bsrbr-eular-moderate-ggamma.csv ttd-bsrbr-eular-moderate-risk.csv
TD_DEPS += ttd-corrona-allpat.csv ttd-corrona-allpat-risk.csv
TD_TARGET := $(DATA)/ttd.rda figs/ttd-corrona-ipdsurv-gengamma.pdf
TD_TARGET += figs/ttd-bsrbr-ipdsurv-comp-gengamma-by-eular.pdf
TD_TARGET += tables/ttd-corrona-ipdsurv-ic.csv tables/ttd-bsrbr-ipdsurv-eular-ic.csv
TD_TARGET += tables/ttd-bsrbr-ipdsurv-eular-adjusted-ic.csv
$(TD_TARGET): $(TD_DEPS)
	Rscript ttd.R
	
# Serious infections
$(DATA)/ttsi.rda:	ttsi.R treatments.csv
	Rscript ttsi.R
	
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
	
# Utility
UTILITY_DEPS := utility.R coef-alava-2013.csv vcov-alava-2013.csv $(DATA)/treatments.rda
UTILITY_TARGET := $(DATA)/pain.rda $(DATA)/utility-mixture.rda $(DATA)/utility-tx-attr.rda
UTILITY_TARGET += $(DATA)/utility-wailoo.rda 
$(UTILITY_TARGET): $(UTILITY_DEPS)
	Rscript utility.R 
	
# Costs
$(DATA)/tx-cost.rda $(DATA)/hosp-cost.rda $(DATA)/mgmt-cost.rda $(DATA)/prod-loss.rda: cost.R cost.xlsx
	Rscript cost.R
	
# Incidence
$(DATA)/incidence-female.rda $(DATA)/incidence-male.rda: incidence.R incidence.csv
	Rscript incidence.R
	