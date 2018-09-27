---
layout: "post"
---
## Variable Composition of Microbial Biomass Required
^1^Petr Capek*, ^1^Katherine Todd-Brown, ^1^Natalie Sadler, ^1^Nancy Hess, ^1^Kirsten Hofmockel

_^1^Pacific Northwest National Laboratory, Environmental Molecular Sciences Laboratory, Richland, WA USA_

## Abstract

## Introduction

## Material and Methods
###Experimental design
Soil microbial culture (see later) was incubated in liquid batch culture microcosms at 24&deg;C in the dark for five days under six different experimental treatments. These treatments include two different organic substrates as a sole source of carbon for microbial culture and three different levels of microcosm structural complexity. The treatments are in full factorial design. The two organic substrates were glucose and cellobiose. The three levels of structural complexity of the system were represented by **_??34 ml??_** incubation vials with liquid medium (denoted as **"BROTH"** further in the text), with liquid medium and mixture of 2.7, 0.1 and 0.1 mm glass beads (5 g of 2.7 mm Biospec glass beads, 2 g of 1.0 mm Biospec glass beads and 1 g of 0.1 mm Next-Advanced glass beads; denoted as **"GLASS"** further in the text) and with liquid medium and glass wool (0.2 g of Pyrex fiber glass wool cutted to ~0.25 cm pieces; denoted as **"WOOL"** further in the text). All experimental treatments were replicated four times **_but this is probably not true since the sampling was destructive. Could you please describe it?_**.

###Soil microbial culture
Soil microbial culture used in the experiment was obtained by imbibement of the **_??prosser soil?? (is it correct)_**, _short characteristic of soil_ ... with the growth medium. One gram of air dried soil was supplemented by the 10 ml of sterile Cellulose Degrader M9 Media (see later) with 6.25 ml of lysogeny broth. The final organic carbon concentration of added medium was ~40 mmol l^-1^ (0.4 mmol g(DW)^-1^). The soil with medium was incubated seven days at 24&deg;C in dark. At the end of incubation, soil-medium suspension was shaken with ~ten 2.7 mm glass beats for ten minutes. The suspension was transferred to sterile flask and soil particles were allowed to settle down. The liquid phase was further transferred to new sterile flask and used as a microbial culture inoculum.

###Microcosm
**_??34 ml??_** empty incubation vials or vials with glass beads and glass wool were filled with 3 ml of sterile medium and 1 ml of microbial culture. The medium was composed of Cellulose Degrader M9 Media and organic substrate (the final concentration 1 g of organic substrate per one liter of medium). Cellulose Degrader M9 Media was prepared by mixing 300 ml of M9 Minimal Media 10X stock solution (59.623 g anhydrous $NaH_{2}PO_{4}$, 29.938 g $KH_{2}PO_{4}$, 4.967 g NaCl, 10.003 g $NH_{4}Cl$, 1.204 g anhydrous $MgSO_{4}$ and 0.140 g $CaCl_{2}\cdot2H_{2}0$ in 1 l of ultra pure water), 10 ml of Hutner's Trace Element Solution (10 g $C_{6}H_{9}NO_{6}$, 7.3 g KOH, 14.45 g $MgSO_{4}$, 3.335 g $CaCl_{2}\cdot2H_{2}0$, 0.00925 g $(NH_{4})_{6}Mo_{7}O_{24}\cdot4H_{2}0$, 0.099 g $FeSO_{4}\cdot7H_{2}0$ and 50 mL of Hutner's Stock Salt Solution in 1 l of ultra pure water), 1 ml of Biotin (1 g $l^{-1}$), 1 ml of Thiamin (1 g $l^{-1}$) and ultra pure water to 1 l. Accounting for the dilution by the inoculum with assumed organic carbon concentration nearly zero, initial organic carbon concentration of the microcosms was 25 mmol $l^{-1}$. Incubation vials were kept covered by Breathe-Easy sealing membrane secured with the aluminium cap at all times the respiration rate was not measured. **_(short decription of the incubator)_**.

###Respiration rate
Microbial respiration rate was measured at the beginning of the experiment and each following day (six times in total). Approximately two hours before the headspace $CO_{2}$ concentration measurement, Breathe-Easy sealing membranes were removed from the incubation vials and replaced by the ruber septa secured with the alluminium cap. Vials headspace was exchanged for the $CO_{2}$-free air and incubated. Headspace was sampled using 5 ml gas tight syringe. 1.5 ml of well mixed headspace air was sampled and directly injected to Li-Cor Li 7000 (LI-COR, Inc., Lincoln, Nebraska, USA). $CO_{2}$ concentration was calculated against the calibartion standard gas with the $CO_{2}$ concentration 2000 ppm. $CO_{2}$ concentration was corrected for the $CO_{2}$ dissolution in liquid media according to Sparling and West (@Sparling1990). Respiration rate was calculated as corrected $CO_{2}$ concentration divided by the time between the headspace atmosphere exchange and the measurement.

###Cellular protein quantification
To quantify microbial biomass in the microcosms, protein content of microbial cells was isolated on the first, second and last day of incubation. The sampling was destructive. The whole volume of microcosm was quantitatively transfered to 15 ml falcon tubes. Cell scraper was used to detach microbial cells from surface of vials or from glass beads/wool. Two ml of phosphate buffer was added and the final solution thoroughly mixed. The final solution was centrifuged at 4,700 g for ten minutes at 4&deg;C and pelet was separated from supernatant. Pellet was resuspended in 0.5 ml and kept deep frozen (-80&deg;C) until the assay was conducted.
Cell pelet was bead beated in ethanol solution for 5 minutes to lyse the cells.Protein concentration was assesed in cell lysate by bicinchoninic acid assay. **_(specification of spectrophotometer)_** was used to measure protein concentration. The concentration was expressed per C basis using the conversion factor 0.45 (@Vrede2004).


###Mathematical description
Three different models (Fig. 1) were used to predict the changes in respiration rate ($R_{H}$) and microbial biomass ($C_{MB}$) in time. The Monod model (Fig. 1a) is a fundamental part of microbial explicit biogeochemical models (@Allison2010; @Wieder2014; @Wieder2013). The Constant Biomass Composition (CBC) model is its adaptation as implemented in MEND model (@Wang2013). The Variable Biomass Composition (VBC) model is an simplified adaptation of Dynamic Energy Budget theory derived by Hanegraaf and Muller (@Hanegraaf2001) for microbial populations. Here we present all models per molar C basis and use the terminology adopted from microbial explicit biogeochemical models for consistency.

####Monod model (Fig. 1a)
In this model, organic substrate ($C_{S}$) is consumed by microbial biomass and transformed to $C_{MB}$ or respired. Microbial biomass is dying at constant rate returning used organic carbon to $C_{S}$ pool:

[1] $~~~~~~~~~~\frac{dC_{MB}}{dt}~=~uptake~\times~CUE~-~k_{MB}~\times~C_{MB}$,

[2] $~~~~~~~~~~\frac{dC_{S}}{dt}~=~-uptake~+~k_{MB}~\times~C_{MB}$.


In eqs. 1 and 2, CUE is the carbon use efficiency (defined as an amount of C incorporated to biomass over the amount of C taken up) and $k_{MB}$ is the death rate constant of microbial biomass decay process. Microbial carbon uptake is defined as hyperbolic function of $C_{S}$ and $C_{MB}$:

[3] $~~~~~~~~~~uptake~=\frac{V_{MAX}~\times~C_{S}~\times~C_{MB}}{K_{M}~+~C_{S}}$,

in which $V_{MAX}$ is maximum velocity constant and $K_{M}$ is affinity constant. Respiration rate is defined as a growth respiration rate ($R_{G}$), which is the constant fraction of uptake:


[4] $~~~~~~~~~~R_{H}~=~R_{G}~=~uptake~\times~(1~-~CUE)$.


####Constant Biomass Composition (CBC) model (Fig. 1b)
In contrast to Monod model, $R_{H}$ consists of two different processes, growth respiration and maintenance respiration ($R_{M}$). $R_{G}$ and $R_{M}$ are respectively defined by following equations:


[5] $~~~~~~~~~~R_{G}~=~(\frac{1}{CUE}~-~1)~\times~\frac{V_{MAX}~\times~C_{S}~\times~C_{MB}}{K_{M}~+~C_{S}}$,

and


[6] $~~~~~~~~~~R_{M}~=~(\frac{1}{CUE}~-~1)~\times~\frac{m_{R}~\times~C_{S}~\times~C_{MB}}{K_{M}~+~C_{S}}$.


In eq. 6, $m_{R}$ is biomass specific maintenance rate constant. Since $R_{G}$ and $R_{M}$ represent a constant fraction of organic C uptake, the equation defining an uptake rate has to take both fluxes into an account and the eq. 3 is rewritten to:


[7] $~~~~~~~~~~uptake~=~\frac{1}{CUE}~\times~(V_{MAX}~+~m_{R}~)\times~\frac{C_{S}~\times~C_{MB}}{K_{M}~+~C_{S}}$.


The overall mass balance equations for pool $C_{MB}$ and $C_{S}$ are respectively:

[8] $~~~~~~~~~~\frac{dC_{MB}}{dt}~=~uptake~-~(R_{G}+R_{M})~-~m_{R}~\times~C_{MB}$,

[9] $~~~~~~~~~~\frac{dC_{S}}{dt}~=~-uptake~+~m_{R}~\times~C_{MB}$.

The last term of eqs. 8 and 9 indicates that the rate of $C_{MB}$ loss due to decay is controled by $m_{R}$.

####Variable Biomass Composition (VBC) model (Fig. 1c)
In contrast to previous models, microbial biomass consits of two functional parts, reserves ($R$) and structures ($S$). Whereas $S$ is functionally similar to $C_{MB}$, $R$ represents a "buffering zone" of $C_{MB}$. $C_{S}$ is assimilated into reserves by a rate proportional to $C_{S}$ and $S$:


[10] $~~~~~~~~~~assimilation~=\frac{V_{MAX}~\times~C_{S}~\times~S}{K_{M}~+~C_{S}}$.


Here we assume that the $C_{S}$ is assimilated into R with efficiency equal to one (@Hanegraaf2001). Organic C is released from reserves by a rate proportional to R. The overall mass balance equation for pool R is thus:

[11] $~~~~~~~~~~\frac{dR}{dt}~=~assimilation~-~f_{0}~\times~R$,

where $f_{0}$ is the constant controlling the rate of release of organic C from reserves. The released organic C is used to maintain structures or to grow (i.e. increase pool $S$). However, maintaining structures have a priority over the growth (@Hanegraaf2001). The growth is realized only when reserves contain enough organic C. When pool $R$ doesn't contain enough organic C to maintain $S$, proportional part of $S$ is lost via respiration. The respective mass balance equations for pools $S$ and $C_{S}$ are defined as:

[12] $~~~~~~~~~~\frac{dS}{dt}~=~max\left\{f_{a}~\times~Y_{S},~0\right\}~+~min\left\{f_{a},~0\right\}$,

and

[13] $~~~~~~~~~~\frac{dC_{S}}{dt}~=~-assimilation$.

In eq. 12, $Y_{S}$ is yield or efficiency of structures construction (functionaly similar to a $CUE$) and $f_{a}$ is the organic C flux available to construct $S$. Depending on the size of the pool $R$ and the amount of organic C needed to maintain $S$, this flux can be positive or negative. If the flux is positive, $S$ increases, while if it is negative, $S$ is lost via respiration (i.e. to cover required maintenance costs). This flux $f_{a}$ is defined as:

[14] $~~~~~~~~~~f_{a}~=~f_{0}~\times~R~-~m_{R}~\times~S$.

$R_{H}$ consists of two different processes, growth and maintenance respiration:

[15] $~~~~~~~~~~R_{H}~=~R_{G}~+~R_{M}~=~~max\left\{f_{a}~\times~(1-Y_{S}),~0\right\}~+~m_{R}~\times~S$.

###Models evaluation
All three models include parameters, whose value can be adjusted to maximize the correspondence between predictions and observations. Since the objective of microbial explicit biogeochmical models is to predict the rate of loss of organic C, we first calibrated the models against the measured respiration rate. Model parameters were adjusted so to minimize the objective function $J$ using the Differential Evolution algorithm (@Mullen2011). The objective function $J$ was defined as:

[16] $~~~~~~~~~~J~=~\sum\limits_{i=1}^{n}{(\frac{P_{i}~-~O_{i}}{\mu})^{2}}$,

where $O_{i}$ and $P_{i}$ stand for observation i and its corresponding value predicted by the model, and $\mu$ is the mean of all observations. Uncertainty of parameters estimates were determined by Constrained Markov Chain Monte Carlo simulation on 5000 iterations (@Soetaert2010). To evaluate the goodness of correspondence between models predictions and observations, log likelihood, Akaike Information Criterion (AIC) and coeficient of determination ($R^{2}$) were calculated. To evaluate the effect of experimental treatments on models prediction capability, model parameters and corresponding goodness of fit were calculated across all treatments or for each substrate (glucose and cellobiose), each level of structural complexity (BROTH, GLASS, WOOL) or each experimetal treatment (combination of the substrate and level of structural complexity) separately.
Second, models were calibrated against measured respiration rate and cellular protein concentration. In this case, the objective function $J$ was defined as a sum of $J$ calculated for respiration rate and cellular protein concentration. The protein concentration was either assumed to be 55% of $C_{MB}$ and 61 and 71% of $R$ and $S$ respectively, or its proportion to $C_{MB}$, $R$ and $S$ was estimated. All analyses were done in statistical program R (@RDevelopmentCoreTeam2014).

## Results

## Discussion

## Figures
![](assets/markdown-img-paste-20180926175020880.png)

$~$
$~$
$~$

![](assets/markdown-img-paste-20180926175259177.png)
##References
