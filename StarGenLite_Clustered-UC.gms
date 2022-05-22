$Title Individual and Clustered UC

$OnText

Developed by

   Diego Tejada (based on model by German Morales-Espana and Andres Ramos)
   Instituto de Investigacion Tecnologica
   Escuela Tecnica Superior de Ingenieria - ICAI
   UNIVERSIDAD PONTIFICIA COMILLAS
   Alberto Aguilera 23
   28015 Madrid, Spain

To Do List

**
$OffText

*-------------------------------------------------------------------------------
*                              MODEL OPTIONS
*-------------------------------------------------------------------------------

* To use a algebraic formulation (formulation without defining the sets)
$OnEmpty OnMulti OffListing

* checking user 1 definition
$if %gams.user1% == "" $log user1 not defined
$if %gams.user1% == "" $stop

* definition of symbol for comments at the end of the line
$EOLCOM //

* optimizer definition
option   lp = cplex ;
option  mip = cplex ;
option rmip = cplex ;

* general options
option optcr    =   1e-4 ;   // tolerance to solve MIP until IntGap < OptcR
option reslim   =  43200 ;   // maximum run time [sec]
option threads  =     -1 ;   // number of cores
option solprint =    off ;   // print the final solution in the .lst file
option limrow   =      0 ;   // maximum number of equations in the .lst file
option limcol   =      0 ;   // maximum number of variables in the .lst file
option bratio   =      1 ;   // bratio = 0 forces GAMS to always try to construct a basis
option savepoint=      2 ;   // save into a gdx file solution (0=no save, 1=only the last one, 2=for each solve)

* profile options
option profile=1, profileTol = 0.01 ;

*-------------------------------------------------------------------------------
*                     MODEL SETS, PARAMETERS, VARIABLES
*-------------------------------------------------------------------------------

* definitions
sets
* periods and scenario sets
   p         "periods e.g.,hours                    "
   pa  (p)   "active   periods                      "
   p1  (p)   "first period of the day               "
   sc        "scenario                              "

*generator sets
   g         "generating         unit               "
   t   (g  ) "individual thermal unit               "
   t1  (g  ) "           thermal units with pMinTU=1"
   ct  (g  ) "clustered  thermal unit               "
   ta  (g  ) "active     thermal unit               "
   tr  (g  ) "thermal units to apply ramping const  "
   tct (g,g) "idividual to clustered thermal unit   "

*Network sets
   i           "bus                                             "
   c           "circuit ID                                      "
   iic (i,i,c) "lines in service in steady state                "
   ii2 (i,i  ) "pair of busses with at least one line connection"
   is  (i    ) "busses without slack bus                        "
   gi  (g  ,i) "generator connected to a bus                    "

   alias (p,pp), (i,ii,j,jj)

parameters

   pTest                   "test" / 0 /
   pTestMinUT(g)           "test min up time clusterizado"

   pExtraClustConst        "extra constraints in the Clustered UC    " / 0 /
   p2ndResUpCost           "cost factor of the up   secondary reserve" / 0.4 /     //0.4
   p2ndResDwCost           "cost factor of the down secondary reserve" / 0.2 /     //0.2
   pDemand        (p     ) "hourly load                    [GW]      "
   pOperReserve   (p     ) "hourly operating reserve       [GW]      "
   pDemand2ndResUp(p     ) "hourly operating reserve up    [GW]      "
   pDemand2ndResDw(p     ) "hourly operating reserve down  [GW]      "
   pIntermGen     (sc,p,i) "stochastic IG generation       [GW]      "
   pScenProb      (  sc  ) "probability of scenarios       [p.u.]    "
   pCommitt       (*,g, p) "commitment  of the unit        [0-1]     "
   pStartup       (  g, p) "startup     in subperiods      [0-1]     "
   pShutDown      (  g, p) "shutdown    in subperiods      [0-1]     "

   pMaxProd       (g     ) "maximum output                 [GW]    "
   pMinProd       (g     ) "minimum output                 [GW]    "
   pIniOut        (g     ) "initial output > min load      [GW]    "
   pIniUC         (g     ) "initial commitment             [0-1]   "
   pIniState      (g     ) "initial state                  [h]     "
   pIniUp         (g     ) "previous up   hours            [h]     "
   pIniDw         (g     ) "previous down hours            [h]     "
   pSUCap         (g     ) "startup  capability            [GW]    "
   pSDCap         (g     ) "shutdown capability            [GW]    "
   pRampUp        (g     ) "ramp up                        [GW/h]  "
   pRampDw        (g     ) "ramp down                      [GW/h]  "
   pMinTU         (g     ) "minimun time down              [h]     "
   pMinTD         (g     ) "minimun time up                [h]     "
   pVariableCost  (g     ) "slope     variable cost        [M�/GWh]"
   pFixCost       (g     ) "               fix cost        [M�/  h]"
   pShutdownCost  (g     ) "shutdown           cost        [M�]    "
   pWeight        (g     ) "cluster weight                 [#]     "
   pENSCost                "energy not served  cost        [M�/GWh]"

   pDurationSD    (g     ) "shutdown duration time         [h]     "
   pDurationSU    (g     ) "startup  duration time         [h]     "
   pTMinSU        (g     ) "min duration time for SU type  [h]     "
   pStartupCost   (g     ) "startup cost                   [M�]    "
   pNumSUTypes    (g     ) "number of active SU types              "
   pWindCurtCost           "cost of wind curtailment       [M�/GWh]"

*Network Parameters
   pSlack                    "slack bus. Bus 1 as default"  / 1 /
   pSBase                    "S Base                                                         [GW]  "
   pNetworkConst             "binary parameter to consider or not the network constraints    [0,1] "
   pBusProfile   (i       )  "load profile per bus                                           [p.u.]"
   pYBUS         (i,ii    )  "Susceptance matrix                                             [p.u.]"
   pYBUSInv      (i,ii    )  "Susceptance matrix inverse Z = B^-1                            [p.u.]"
   pPmax         (i,ii,c  )  "maximum Power Transfer of line                                 [GW]  "
   pISF          (i,ii,c,i)  "injection Shift Factors (for node i to line i->ii circuitID c) [p.u.]"
   pXline        (i,ii,c  )  "Reactance  X of line                                           [p.u.]"
   pRline        (i,ii,c  )  "Resistance R of line                                           [p.u.]"
   pZline        (i,ii,c  )  "Impedance  Z of line                                           [p.u.]"

*Output Parameters
   pSummary      (* ,*       ) "Results summary                                "
   pProduct      (*,sc,g,p   ) "output      of the unit                 [MW]   "
   pIG           (  sc,p,i   ) "intermittent generation                 [MW]   "
   pENS          (*,sc,p,i   ) "energy not served per node              [MW]   "
   pSRMC         (  sc,p,i   ) "short run marginal cost                 [�/MWh]"
   p2ndResUP     (*, g,p     ) "up   secondary reserve                  [MW]   "
   p2ndResDW     (*, g,p     ) "down secondary reserve                  [MW]   "
   p2ResUpMC     (     p     ) "marginal cost of up   secondary reserve [�/MW] "
   p2ResDwMC     (     p     ) "marginal cost of down secondary reserve [�/MW] "
   pPrintCirPF   (sc,  i,ii,c) "Lines to be printed in steady state            "
   pCirPF        (sc,p,i,ii,c) "Power flow in steady state              [MW]   "

variables
   vTotalVCost             "total system variable cost      [M�]"
   vCirPF    (sc,p,i,ii,c) "power flow through a line       [GW]"

integer  variables
   vCommit       (   p,g)  "commitment of the unit         [0-1]"
   vStartUp      (   p,g)  "startup    of the unit         [0-1]"
   vShutDown     (   p,g)  "shutdown   of the unit         [0-1]"

positive variables
   vProduct      (sc,p,g)  "output of the unit              [GW]"
   vProduct1     (sc,p,g)  "output of the unit > minload    [GW]"
   v2ndResUP     (   p,g)  "2nd res. up   allocation        [GW]"
   v2ndResDW     (   p,g)  "2nd res. down allocation        [GW]"
   vIG           (sc,p,i)  "intermittent generation per bus [GW]"
   vENS          (sc,p,i)  "energy not served per node      [GW]"

*-------------------------------------------------------------------------------
*                                EQUATIONS
*-------------------------------------------------------------------------------

equations
   eTotalVCost              "tot system variable cost                   [M�]"
   eUCStrShut  (   p,g    ) "relation among committment startup and shutdown"
   eMinOutput  (sc,p,g    ) "min output of a committed unit             [GW]"
   eMinTUp     (sc,p,g    ) "minimun up time   (    committed)              "
   eMinTDw     (sc,p,g    ) "minimun down time (not committed)              "
   eRampUp     (sc,p,g    ) "bound on ramp up                           [GW]"
   eRampDw     (sc,p,g    ) "bound on ramp down                         [GW]"
   eTotOut     (sc,p,g    ) "tot output of a        committed unit      [GW]"
   eMaxOut     (sc,p,g    ) "max output of a        committed unit      [GW]"
   eMaxOut2    (sc,p,g    ) "max output of a        committed unit      [GW]"
   eMaxOut3    (sc,p,g    ) "max output of a        committed unit      [GW]"
*Extra Clustered UC constraints
   eCUC        (   p,g    ) "relation among individual and clustered units  "
   eMaxOut4    (sc,p,g    ) "max output of a        committed unit      [GW]"
   eTotOut2    (sc,p,g    ) "tot output of a clustered unit             [GW]"
   eTotResUp   (   p,g    ) "tot reserve up   of a clustered unit       [GW]"
   eTotResDw   (   p,g    ) "tot reserve down of a clustered unit       [GW]"
   eNonSym1    (   p,g,g  ) "constraint to avoid symetry                    "
   eNonSym2    (   p,g,g  ) "constraint to avoid symetry                    "
   eNonSym3    (   p,g,g  ) "constraint to avoid symetry                    "
*Reserve constraints
   e2ReserveUp (   p      ) "operating reserve upwards                  [GW]"
   e2ReserveDw (   p      ) "operating reserve downwards                [GW]"
*Network constraints
   eBalanceBus (sc,p,i    ) "load generation balance per bus            [GW]"
   eBalance    (sc,p      ) "load generation balance                    [GW]"
   eCirPF      (sc,p,i,i,c) "power flow through a line                  [GW]"
;

*--------------------------- Objective Function --------------------------------

eTotalVCost ..
   vTotalVCost =e=
   + sum[(sc,pa(p),i    ) $[pScenProb(sc)], pENSCost                       * vENS     (sc,p,i)]
   + sum[(sc,pa(p),ta(g)) $[pScenProb(sc)], pVariableCost(g)               *[vProduct1(sc,p,g)
                                                               +pMinProd(g)* vCommit  (   p,g)]]
   + sum[(   pa(p),ta(g))                 , pFixCost     (g)               * vCommit  (   p,g)]
   + sum[(   pa(p),ta(g))                 , pStartupCost (g)               * vStartUp (   p,g)]
   + sum[(   pa(p),ta(g))                 , pShutdownCost(g)               * vShutDown(   p,g)]
   + sum[(   pa(p),ta(g))                 , pFixCost     (g)*p2ndResUpCost * v2ndResUP(   p,g)]
   + sum[(   pa(p),ta(g))                 , pFixCost     (g)*p2ndResDwCost * v2ndResDW(   p,g)]
   + sum[(sc,pa(p),i    ) $[pScenProb(sc)], pWindCurtCost *[pIntermGen(sc,p,i)-vIG    (sc,p,i)]]
;

*---------------------------- UC Formulation -----------------------------------

eUCStrShut    (pa(p),ta(g)) ..
   vCommit (p,g) - vCommit (p-1,g) - pIniUC(g) $p1(p) =e= vStartUp(p,g) - vShutDown(p,g) ;

eMinOutput    (sc,pa(p),ta(g)) $[pScenProb(sc)] ..
   vProduct1(sc,p,g) - v2ndResDW(p,g) =g= 0 ;

eMinTUp       (sc,pa(p),ta(g)) $[pScenProb(sc) and ord(p)>=pMinTU(g)] ..
   sum[pp $[ord(pp)>=ord(p)+1-pMinTU(g) and ord(pp)<=ord(p) and pDemand(pp)], vStartUp (pp,g)] =l=              vCommit(p,g) ;

eMinTDw       (sc,pa(p),ta(g)) $[pScenProb(sc) and ord(p)>=pMinTD(g)] ..
   sum[pp $[ord(pp)>=ord(p)+1-pMinTD(g) and ord(pp)<=ord(p) and pDemand(pp)], vShutDown(pp,g)] =l= pWeight(g) - vCommit(p,g) ;

eRampUp       (sc,pa(p),tr(g)) $[pScenProb(sc)] ..
   + vProduct1(sc,p  ,g)
   - vProduct1(sc,p-1,g)
   + v2ndResUP(   p  ,g)
   - [pIniOut(g)-pMinProd(g)*pIniUC(g)] $[p1(p) and pIniUC(g)]
  =l=
   + pRampUp (g) * vCommit(p,g)
;

eRampDw       (sc,pa(p),tr(g)) $[pScenProb(sc)] ..
   + vProduct1(sc,p  ,g)
   - vProduct1(sc,p-1,g)
   - v2ndResDW(   p  ,g)
   - [pIniOut(g)-pMinProd(g)*pIniUC(g)] $[p1(p) and pIniUC(g)]
  =g=
   - pRampDw (g) *[vCommit(p-1,g) + pIniUC(g)$p1(p)]
;

eTotOut       (sc,pa(p),ta(g)) $[pScenProb(sc)] ..
   + vProduct (sc,p  ,g)
  =e=
   + pMinProd(g)* vCommit  (   p,g)
   +              vProduct1(sc,p,g)
;

eMaxOut       (sc,pa(p),ta(g)) $[pScenProb(sc) and not[t1(g)]] ..
   + vProduct1(sc,p  ,g)
   + v2ndResUP(   p  ,g)
  =l=
   + vCommit  (   p  ,g)*   [pMaxProd(g)-pMinProd(g)  ]
   - vStartUp (   p  ,g)*   [pMaxProd(g)-pSUCap  (g)  ]
   - vShutDown(   p+1,g)*   [pMaxProd(g)-pSDCap  (g)  ]
;

eMaxOut2      (sc,pa(p),ta(g)) $[pScenProb(sc) and     t1(g)] ..
   + vProduct1(sc,p  ,g)
   + v2ndResUP(   p  ,g)
  =l=
   + vCommit  (   p  ,g)*   [pMaxProd(g)-pMinProd(g)  ]
   - vStartUp (   p  ,g)*   [pMaxProd(g)-pSUCap  (g)  ]
   - vShutDown(   p+1,g)*max[pSUCap  (g)-pSDCap  (g),0]
;

eMaxOut3      (sc,pa(p),ta(g)) $[pScenProb(sc) and     t1(g)] ..
   + vProduct1(sc,p  ,g)
   + v2ndResUP(   p  ,g)
  =l=
   + vCommit  (   p  ,g)*   [pMaxProd(g)-pMinProd(g)  ]
   - vShutDown(   p+1,g)*   [pMaxProd(g)-pSDCap  (g)  ]
   - vStartUp (   p  ,g)*max[pSDCap  (g)-pSUCap  (g),0]
;

*------------------- Reserve requirements Constraints --------------------------

e2ReserveUp(pa(p)) $[pDemand2ndResUp(p)].. sum[ta(g), v2ndResUP(p,g)] =g= pDemand2ndResUp(p) ;
e2ReserveDw(pa(p)) $[pDemand2ndResDw(p)].. sum[ta(g), v2ndResDW(p,g)] =g= pDemand2ndResDw(p) ;

*------------------------- Network Constraints ---------------------------------

eBalanceBus(sc,pa(p),i) $[pScenProb(sc) and pNetworkConst=1] ..
    + sum[gi(ta(g),i),vProduct(sc,p,g     )]
    +                 vIG     (sc,p,i     )
    +                 vENS    (sc,p,i     )
    + sum[iic(ii,i,c),vCirPF  (sc,p,ii,i,c)]
    - sum[iic(i,ii,c),vCirPF  (sc,p,i,ii,c)]
   =e=
    + pDemand(p)*pBusProfile(i)
;

eBalance   (sc,pa(p)  ) $[pScenProb(sc) and pNetworkConst=0] ..
    + sum[ta(g), vProduct(sc,p,g)]
    + sum[i    , vIG     (sc,p,i)
    +            vENS    (sc,p,i)]
   =e=
    + sum[i    , pDemand(p)*pBusProfile(i)]
;

eCirPF      (sc,pa(p),iic(i,ii,c)) $[pScenProb(sc) and pNetworkConst=1] ..
   + vCirPF (sc,p,iic)
  =e=
   + sum[is, pISF(iic,is) *
           [+ sum[gi(ta(g),is),vProduct(sc,p,g )]
            +                  vIG     (sc,p,is)
            +                  vENS    (sc,p,is)
            -                  pDemand (   p   ) * pBusProfile(is)
           ]
        ]
;

*--- Extra Constraints to include individual information in the Clustered UC ---

eCUC    (   pa(p),ta(g)) $[pExtraClustConst=1].. vCommit(p,g) =e= sum[tct(t,g),vCommit(p,t)] ;

eMaxOut4(sc,pa(p),t    ) $[pExtraClustConst=1 and pScenProb(sc)]..
   + vProduct1(sc,p  ,t)
   + v2ndResUP(   p  ,t)
  =l=
   + vCommit  (   p  ,t)* [pMaxProd(t)-pMinProd(t)]
;

eTotOut2 (sc,pa(p),  ta(g)) $[pExtraClustConst=1 and pScenProb(sc)]..
   + vProduct1(sc,p ,g ) =e= sum[tct(t,g),vProduct1(sc,p,t)] ;

eTotResUp(   pa(p),  ta(g)) $[pExtraClustConst=1]..
   + v2ndResUP(   p ,g ) =e= sum[tct(t,g),v2ndResUP(   p,t)] ;

eTotResDw(   pa(p),  ta(g)) $[pExtraClustConst=1]..
   + v2ndResDw(   p ,g ) =e= sum[tct(t,g),v2ndResDw(   p,t)] ;

*eNonSym1(   pa(p),t,ta(g)) $[pExtraClustConst=1 and tct(t  ,g) and pMinTU(t)<=pTestMinUT(t)]..
eNonSym1(   pa(p),t,ta(g)) $[pExtraClustConst=1 and tct(t  ,g)]..
   + sum[tct(t+1,g), vCommit(p,t+1)] =l= sum[tct(t,g), vCommit(p,t)] ;

eNonSym2(   pa(p),t,ta(g)) $[pExtraClustConst=2 and tct(t-1,g) and pMinTU(t)<=pTestMinUT(t)]..
   + sum[tct(t-1,g), vCommit(p,t-1)] =l= 1+sum[tct(t,g), vCommit(p,t)]-sum[tct(t,g), vCommit(p-1,t) + pIniUC(t) $p1(p)];

eNonSym3(   pa(p),t,ta(g)) $[pExtraClustConst=2 and tct(t-1,g) and pMinTU(t)> pTestMinUT(t)]..
   + sum[tct(t-1,g), vCommit(p,t-1)] =g=   sum[tct(t,g), vCommit(p,t)]-sum[tct(t,g), vCommit(p-1,t) + pIniUC(t) $p1(p)];

*-------------------------------------------------------------------------------
*                                MODELS
*-------------------------------------------------------------------------------

model SDUC  / all / ;
 SDUC.holdfixed = 1 ;  SDUC.optfile = 1 ;  SDUC.trylinear = 1 ;

model RSDUC / all / ;
RSDUC.holdfixed = 1 ; RSDUC.optfile = 1 ;

*-------------------------------------------------------------------------------
*                             DATA FROM EXCEL FILE
*-------------------------------------------------------------------------------

* read input data from Excel and include into the model
file TMP / tmp.txt /
$onecho  > tmp.txt
   i="%gams.user1%.xlsx"
   r1=indices
   o1=indices
   r2=param
   o2=param
   r3=demand
   o3=demand
   r4=oprres
   o4=oprres
   r5=IGgen
   o5=IGgen
   r6=thermalgen
   o6=thermalgen
   r7=demandbus
   o7=demandbus
   r8=network
   o8=network
   r9=clusteredgen
   o9=clusteredgen
   r10=gentocluster
   o10=gentocluster
$offecho
$call xls2gms m @"tmp.txt"

sets
$include  indices
$include  gentocluster
;
$include  param
;

parameter pBusProfile(i)        load profile per bus   [p.u.] /
$include  demandbus
                                                              /

table     tDemand(p,*)          hourly    load           [MW]
$include  demand

table     tDemand2ndRes(p,*)    hourly operating reserve [MW]
$include  oprres

table     tIntermGen(p,i,*,sc)  stochastic IG generation [MW]
$include  IGgen

table     tThermalGen(g,*)      individual thermal generation data
$include  thermalgen

table     tClusteredGen(g,*)    clustered  thermal generation data
$include  clusteredgen

table     tNetwork(i,ii,c,*)    Network parameters
$include  network

execute 'del tmp.txt indices param demand oprres demandbus IGgen ' ;
execute 'del tmp.txt thermalgen network clusteredgen gentocluster' ;

*-------------------------------------------------------------------------------
*                   OPTIONS DEFINITION FOR SOLVERS
*-------------------------------------------------------------------------------

FILE     GOPT / gurobi.opt /
*PUT      GOPT / 'IIS 1'    / 'lazyconstraints 1' / 'eCUC.lazy 3' /
PUT      GOPT / 'IIS 1'    /
PUTCLOSE GOPT ;

FILE     COPT / cplex.opt  /
PUT      COPT / 'IIS yes'  /
PUTCLOSE COPT ;

*-------------------------------------------------------------------------------
*                      Injection Sensitivity Factors
*-------------------------------------------------------------------------------

*Busses without the slack bus
is(i) $[ORD(i) <> pSlack] = YES;

if(pNetworkConst=1,
*Line conecttions
   iic (i,ii,c) $ tNetwork(i,ii,c,'InService') = YES                ;
   ii2 (i,ii  )                                = SUM[c,iic(i,ii,c)] ;

*Obtaining the parameters
   pXline (i,ii,c) = tNetwork(i,ii,c,'X');

*Obtaining susceptance matrix with multiple circuits defined as i->j and j->i
   pYBUS(ii2(i,ii)) = -SUM[c $ (iic(i,ii,c) AND pXline (i,ii,c)),1/pXline(i,ii,c)]
                      -SUM[c $ (iic(ii,i,c) AND pXline (ii,i,c)),1/pXline(ii,i,c)];

*Creating the symmetric matrix
   pYBUS(ii,i) $ii2(i,ii) = pYBUS(i,ii);

*Creating the diagonal for the B matriz
   pYBUS(ii,ii) = SUM(i,-pYBUS(i,ii));

*Eliminaing the Slack/reference node
   pYBUS(i,ii) $[NOT(is(i) )] = 0;
   pYBUS(i,ii) $[NOT(is(ii))] = 0;

*Creating gdx files for inverse calculation
   EXECUTE_UNLOAD 'gdxfrominverse.gdx' is ;
   EXECUTE_UNLOAD 'gdxforinverse.gdx ' is ;

*Obtaining the inverse of pYBUS and saving it in pYBUSInv
   EXECUTE_UNLOAD 'gdxforinverse.gdx' is,pYBUS;
   EXECUTE 'invert gdxforinverse.gdx  is pYBUS gdxfrominverse.gdx pYBUSInv';
   EXECUTE_LOAD  'gdxfrominverse.gdx' ,  pYBUSInv;

*Obtaining the Injection Shift Factors
   pISF(i,ii,c,is) $ [pXline(i,ii,c) AND iic(i,ii,c)]
                   = (pYBUSInv(i,is) - pYBUSInv(ii,is))/pXline(i,ii,c);
*Eliminating the pISF lower than 1E-9 (contribution lower than 1W)
   pISF(i,ii,c,is) $ [ABS[pISF(i,ii,c,is)] < 1E-9] = 0;
);

*-------------------------------------------------------------------------------
*                        SET VALUE VARIABLES AND DYNAMIC SETS
*-------------------------------------------------------------------------------

* scaling of demand and renewable energy to GW
pDemand        (p) = tDemand      (p,'Energy'   ) * 1e-3 ;
pDemand2ndResUp(p) = tDemand2ndRes(p,'Energy_up') * 1e-3 ;
pDemand2ndResDw(p) = tDemand2ndRes(p,'Energy_dw') * 1e-3 ;
pIntermGen(sc,p,i) = tIntermGen   (p,i,'En',sc  ) * 1e-3 ;

* active periods with demand
pa(p) $[pDemand(p)] = yes ;

* determine the first hour of the day
p1(p) $[ord(p) = 1] = yes ;

* assignment of thermal units and active thermal units
*t (g) $[tThermalGen  (g,'FuelCost') * tThermalGen  (g,'MaxProd')] = yes ;
*ct(g) $[tClusteredGen(g,'FuelCost') * tClusteredGen(g,'MaxProd')] = yes ;
ta(g) $[t(g)] = yes ;
tr(g) $[t(g)] = yes ;

*determine the connection between the generator and the bus
gi( t,i) $[tThermalGen  ( t,'Bus') and [tThermalGen  ( t,'Bus')=ord(i)]] = yes ;
gi(ct,i) $[tClusteredGen(ct,'Bus') and [tClusteredGen(ct,'Bus')=ord(i)]] = yes ;


* scaling of parameters to GW and M�
pENSCost          = pENSCost                             * 1e-3 ;
pWindCurtCost     = pWindCurtCost                        * 1e-3 ;
pMaxProd     (t ) = tThermalGen   (t ,'MaxProd'        ) * 1e-3 ;
pMinProd     (t ) = tThermalGen   (t ,'MinProd'        ) * 1e-3 ;
pIniOut      (t ) = tThermalGen   (t ,'IniProd'        ) * 1e-3 ;
pSUCap       (t ) = tThermalGen   (t ,'SUcap'          ) * 1e-3 ;
pSDCap       (t ) = tThermalGen   (t ,'SDcap'          ) * 1e-3 ;
pRampUp      (t ) = tThermalGen   (t ,'RampUp'         ) * 1e-3 ;
pRampDw      (t ) = tThermalGen   (t ,'RampDw'         ) * 1e-3 ;
pIniState    (t ) = tThermalGen   (t ,'IniState'       ) ;
pMinTU       (t ) = tThermalGen   (t ,'MinTU'          ) ;
pMinTD       (t ) = tThermalGen   (t ,'MinTD'          ) ;
pVariableCost(t ) = tThermalGen   (t ,'OMVarCost'      ) * 1e-3 +
                    tThermalGen   (t ,'SlopeVarCost'   ) * 1e-3 * tThermalGen  (t ,'FuelCost') ;
pDurationSD  (t ) = tThermalGen   (t ,'SDduration'     ) ;
pFixCost     (t ) = tThermalGen   (t ,'InterVarCost'   ) * 1e-6 * tThermalGen  (t ,'FuelCost') ;
pShutdownCost(t ) = tThermalGen   (t ,'ShutdownCost'   ) * 1e-6 * tThermalGen  (t ,'FuelCost') ;
pDurationSU  (t ) = tThermalGen   (t ,'SUduration1'    ) ;
pTMinSU      (t ) = tThermalGen   (t ,'DownTtimeforSU1') ;
pStartupCost (t ) = tThermalGen   (t ,'SUcost1'        ) * 1e-6 * tThermalGen  (t ,'FuelCost') ;

pMaxProd     (ct) = tClusteredGen (ct,'MaxProd'        ) * 1e-3 ;
pMinProd     (ct) = tClusteredGen (ct,'MinProd'        ) * 1e-3 ;
pSUCap       (ct) = tClusteredGen (ct,'SUcap'          ) * 1e-3 ;
pSDCap       (ct) = tClusteredGen (ct,'SDcap'          ) * 1e-3 ;
pRampUp      (ct) = tClusteredGen (ct,'RampUp'         ) * 1e-3 ;
pRampDw      (ct) = tClusteredGen (ct,'RampDw'         ) * 1e-3 ;
pIniState    (ct) = tClusteredGen (ct,'IniState'       ) ;
pMinTU       (ct) = tClusteredGen (ct,'MinTU'          ) ;
pMinTD       (ct) = tClusteredGen (ct,'MinTD'          ) ;
pVariableCost(ct) = tClusteredGen (ct,'OMVarCost'      ) * 1e-3 +
                    tClusteredGen (ct,'SlopeVarCost'   ) * 1e-3 * tClusteredGen(ct,'FuelCost') ;
pDurationSD  (ct) = tClusteredGen (ct,'SDduration'     ) ;
pFixCost     (ct) = tClusteredGen (ct,'InterVarCost'   ) * 1e-6 * tClusteredGen(ct,'FuelCost') ;
pShutdownCost(ct) = tClusteredGen (ct,'ShutdownCost'   ) * 1e-6 * tClusteredGen(ct,'FuelCost') ;
pDurationSU  (ct) = tClusteredGen (ct,'SUduration1'    ) ;
pTMinSU      (ct) = tClusteredGen (ct,'DownTtimeforSU1') ;
pStartupCost (ct) = tClusteredGen (ct,'SUcost1'        ) * 1e-6 * tClusteredGen(ct,'FuelCost') ;

*Network parameters
pSBase       = pSBase                         * 1e-3  ;
pPmax  (iic) = tNetwork(iic,'Pmax')           * 1e-3  ;
pRline (iic) = tNetwork(iic,'R'   )                   ;
pXline (iic) = tNetwork(iic,'X'   )                   ;
pZline (iic) = SQRT[SQR[pRline(iic)]+SQR[pXline(iic)]];

* Weight for individuals is equal to 1, while for clustered depends on the amount of units in the cluster
pWeight(t ) = 1 ;
pWeight(ct) = sum[tct(t,ct),1] ;

* if the initial output of the unit is above its minimum load then the unit is committed, otherwise it is not committed
pIniUC ( t) =                1          $[pIniOut(t) >= pMinProd(t)]  ;
pIniUC (ct) = sum[tct(t,ct), 1          $[pIniOut(t) >= pMinProd(t)]] ;
pIniOut(ct) = sum[tct(t,ct), pIniOut(t) $[pIniOut(t) >= pMinProd(t)]] ;

* if the minimun up or down times are 0, they are changed to 1
pMinTU(g) $[pMinTU(g)=0] = 1 ;
pMinTD(g) $[pMinTD(g)=0] = 1 ;

*Initial up/dw/Start-Up Conditions
pIniUp(g) $[pIniState(g)>=0] = max[0,[pMinTU(g)-    pIniState(g) ]*            pIniUC(g) ];
pIniDw(g) $[pIniState(g)< 0] = max[0,[pMinTD(g)-abs[pIniState(g)]]*[pWeight(g)-pIniUC(g)]];

* definition of units with pMinTU=1
t1(g) $[pMinTU(g)=1] = yes ;

* ensuring that input data are within bounds
pSUCap (g) = min[pMaxProd(g)            ,max[pMinProd(g),pSUCap (g)]] ;
pSDCap (g) = min[pMaxProd(g)            ,max[pMinProd(g),pSDCap (g)]] ;
pRampUp(g) = min[pMaxProd(g)-pMinProd(g),                pRampUp(g) ] ;
pRampDw(g) = min[pMaxProd(g)-pMinProd(g),                pRampDw(g) ] ;

* extra ramps fix cost
pStartupCost (g) = pStartupCost (g) + pDurationSU(g)*pFixCost(g) + pDurationSU(g)*pMinProd(g)*pVariableCost(g)/2;
pShutdownCost(g) = pShutdownCost(g) + pDurationSD(g)*pFixCost(g) + pDurationSD(g)*pMinProd(g)*pVariableCost(g)/2;

* test X=(Pmax+UR-SU)/(UR+DR)
pTestMinUT(t) = [pMaxProd(t)+pRampUp(t)-pSUCap(t)]/[pRampUp(t)+pRampDw(t)];

*display pTestMinUT;

*-------------------------------------------------------------------------------
*                        SET BOUNDS FOR VARIABLES
*-------------------------------------------------------------------------------

* for individual thermal units the upper bound is equal to 1 (binary)
vCommit.up  (pa(p), t) = 1 ;
vStartUp.up (pa(p), t) = 1 ;
vShutDown.up(pa(p), t) = 1 ;

* upper bound on clustered generator depending on the number of individual generators in the cluster
vCommit.up  (pa(p),ct) = pWeight(ct) ;
vStartUp.up (pa(p),ct) = pWeight(ct) ;
vShutDown.up(pa(p),ct) = pWeight(ct) ;

* bounds on variables
vProduct.up (sc,pa(p),g) $[pScenProb(sc)] = pMaxProd    (g     ) * pWeight(g) ;
vProduct1.up(sc,pa(p),g) $[pScenProb(sc)] =[pMaxProd    (g     ) -
                                            pMinProd    (g     )]* pWeight(g) ;
vIG.up      (sc,pa(p),i) $[pScenProb(sc)] = pIntermGen  (sc,p,i) ;
vENS.up     (sc,pa(p),i) $[pScenProb(sc)] = pDemand     (p     ) *
                                            pBusProfile (i     ) ;

* bound for ramps
v2ndResUP.up(pa(p),g) = pWeight(g) * min[pRampUp(g)+pRampDw(g),pMaxProd(g)-pMinProd(g)];
v2ndResDW.up(pa(p),g) = pWeight(g) * min[pRampUp(g)+pRampDw(g),pMaxProd(g)-pMinProd(g)];

v2ndResUP.fx(pa(p),g) $[pDemand2ndResUp(p)=0] = 0 ;
v2ndResDW.fx(pa(p),g) $[pDemand2ndResDw(p)=0] = 0 ;

* bound for initial up/dw/Start-Up Conditions
vCommit.fx(pa(p),g) $[ord(p)<=[pIniUp(g) + pIniDw(g)]] = pIniUC(g);

* bounds for circuit power flow in steady state
vCirPF.up (sc,pa(p),iic  )  $[pScenProb(sc)  ] =  pPmax  (iic);
vCirPF.lo (sc,pa(p),iic  )  $[pScenProb(sc)  ] = -pPmax  (iic);

* if network is not used then we only need one vENS
vENS.fx   (sc,pa(p),is(i))  $[pNetworkConst=0] = 0 ;

*-------------------------------------------------------------------------------
*                           SOLVE INDIVIDUAL UC
*-------------------------------------------------------------------------------
*$ontext
solve RSDUC using RMIP minimizing vTotalVCost ;
solve  SDUC using  MIP minimizing vTotalVCost ;

pSummary('Final Obj. Function           ','Individual-UC') = SDUC.Objval   + eps ;
pSummary('CPU Time  Model generation [s]','Individual-UC') = SDUC.resGen   + eps ;
pSummary('CPU Time  Model solution   [s]','Individual-UC') = SDUC.resUsd   + eps ;
pSummary('Number of nodes used in MIP   ','Individual-UC') = SDUC.nodUsd   + eps ;
pSummary('Number of variables           ','Individual-UC') = SDUC.numVar   + eps ;
pSummary('Number of discrete variables  ','Individual-UC') = SDUC.numDVar  + eps ;
pSummary('Number of equations           ','Individual-UC') = SDUC.numEqu   + eps ;
pSummary('Number of nonzero elements    ','Individual-UC') = SDUC.numNZ    + eps ;
pSummary('Integrality gap               ','Individual-UC') =abs[SDUC.Objval-RSDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Final relative gap            ','Individual-UC') =abs[SDUC.Objest- SDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Network constraints considered','Individual-UC') = pNetworkConst + eps ;
pSummary('Fix          Cost'             ,'Individual-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)                 * vCommit.l  (p,g   )] + eps ;
pSummary('Startup      Cost'             ,'Individual-UC') = sum[(p,ta(g)   )                 , pStartupCost (g)                 * vStartUp.l (p,g   )] + eps ;
pSummary('ShutDown     Cost'             ,'Individual-UC') = sum[(p,ta(g)   )                 , pShutdownCost(g)                 * vShutDown.l(p,g   )] + eps ;
pSummary('ReserveUp    Cost'             ,'Individual-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResUpCost   * v2ndResUP.l(p,g   )] + eps ;
pSummary('ReserveDw    Cost'             ,'Individual-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResDwCost   * v2ndResDW.l(p,g   )] + eps ;
pSummary('ENS          Cost'             ,'Individual-UC') = sum[(sc,p,i    ) $[pScenProb(sc)], pENSCost                         * vENS.l     (sc,p,i)] + eps ;
pSummary('Variable     Cost'             ,'Individual-UC') = sum[(sc,p,ta(g)) $[pScenProb(sc)], pVariableCost(g)*[vProduct1.l(sc,p,g) + pMinProd(g)* vCommit.l(   p,g)]] + eps ;
pSummary('Curtailment  Cost'             ,'Individual-UC') = sum[(sc,p,i    ) $[pScenProb(sc)], pWindCurtCost   *[pIntermGen (sc,p,i)               -vIG.l    (sc,p,i)]] + eps ;

* scaling of results

* individual results
pCommitt ('Individual-UC',    t,pa(p))                     = vCommit.l    (   p,t)     + eps ;
pProduct ('Individual-UC',sc, t,pa(p))$[pScenProb(sc)    ] = vProduct.l   (sc,p,t)*1e3 + eps ;
* clustered results
pCommitt ('Individual-UC',   ct,pa(p))                     = sum[tct(t,ct), vCommit.l (   p,t)    ] + eps ;
pProduct ('Individual-UC',sc,ct,pa(p))$[pScenProb(sc)    ] = sum[tct(t,ct), vProduct.l(sc,p,t)*1e3] + eps ;

pENS     ('Individual-UC',sc,pa(p),i) $[pScenProb(sc)    ] = vENS.l       (sc,p,i)*1e3 + eps ;
pSRMC    (sc,pa(p),i) $[pScenProb(sc) and pNetworkConst=1] = eBalanceBus.m(sc,p,i)*1e3 + eps ;
pSRMC    (sc,pa(p),i) $[pScenProb(sc) and pNetworkConst=0] = eBalance.m   (sc,p  )*1e3 + eps ;
p2ResUpMC(   pa(p)  )                                      = e2ReserveUp.m(   p  )*1e3 + eps ;
p2ResDwMC(   pa(p)  )                                      = e2ReserveDw.m(   p  )*1e3 + eps ;
p2ndResUP('Individual-UC',ct,pa(p)  )                      = sum[tct(t,ct), v2ndResUP.l(p,t)*1e3] + eps ;
p2ndResDW('Individual-UC',ct,pa(p)  )                      = sum[tct(t,ct), v2ndResDW.l(p,t)*1e3] + eps ;
p2ndResUP('Individual-UC', t,pa(p)  )                      =                v2ndResUP.l(p,t)*1e3  + eps ;
p2ndResDW('Individual-UC', t,pa(p)  )                      =                v2ndResDW.l(p,t)*1e3  + eps ;
pIG      (sc,pa(p),i) $[pScenProb(sc) and pIntermGen(sc,p,i)]= vIG.l      (sc,p,i)*1e3 + eps ;

*Parameters for only print the lines that are congested in steady state
pPrintCirPF(sc      ,iic) $[sum(p,abs(vCirPF.m (sc,p,iic))>0)    ] = 1 ;
pCirPF     (sc,pa(p),iic) $[pScenProb(sc) and pPrintCirPF(sc,iic)] = vCirPF.l(sc,p,iic)*1e3 + eps ;

* data output to xls file
put TMP put 'par=pCommitt  rdim=2 rng=UC!a1'      / 'par=pProduct  rdim=3 rng=Output!a1'  / 'par=pIG       rdim=2 rng=IG!a1'         /
put TMP put 'par=pSRMC     rdim=2 rng=SRMC!a1'    / 'par=p2ndResUP rdim=2 rng=SecResUP!a1'/ 'par=p2ndResDW rdim=2 rng=SecResDW!a1'   / 'par=p2ResUpMC rdim=1 rng=SecResMC!a1' /
put TMP put 'par=p2ResDwMC rdim=1 rng=SecResMC!c1'/ 'par=pCirPF    rdim=2 rng=CirPF!a1'   /
put TMP put 'par=pENS      rdim=3 rng=ENS!a1'     / 'par=pSummary  rdim=1 rng=Summary!a1' /
putclose
execute_unload   'tmp.gdx' pProduct pCommitt pIG pSRMC p2ndResUP p2ndResDW p2ResUpMC p2ResDwMC pCirPF pENS pSummary
execute          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
execute          'del        tmp.gdx                             tmp.txt' ;

*-------------------------------------------------------------------------------
*                           CLUSTERED UC
*-------------------------------------------------------------------------------

* reset and new definition of active thermal units
ta(g)          = no  ;
tr(g)          = no  ;
ta(g) $[ct(g)] = yes ;
tr(g) $[ct(g)] = yes ;

solve RSDUC using RMIP minimizing vTotalVCost ;
solve  SDUC using  MIP minimizing vTotalVCost ;

pSummary('Final Obj. Function           ','Clustered-UC') = SDUC.Objval   + eps ;
pSummary('CPU Time  Model generation [s]','Clustered-UC') = SDUC.resGen   + eps ;
pSummary('CPU Time  Model solution   [s]','Clustered-UC') = SDUC.resUsd   + eps ;
pSummary('Number of nodes used in MIP   ','Clustered-UC') = SDUC.nodUsd   + eps ;
pSummary('Number of variables           ','Clustered-UC') = SDUC.numVar   + eps ;
pSummary('Number of discrete variables  ','Clustered-UC') = SDUC.numDVar  + eps ;
pSummary('Number of equations           ','Clustered-UC') = SDUC.numEqu   + eps ;
pSummary('Number of nonzero elements    ','Clustered-UC') = SDUC.numNZ    + eps ;
pSummary('Integrality gap               ','Clustered-UC') =abs[SDUC.Objval-RSDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Final relative gap            ','Clustered-UC') =abs[SDUC.Objest- SDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Network constraints considered','Clustered-UC') = pNetworkConst + eps ;
pSummary('Fix          Cost'             ,'Clustered-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)                 * vCommit.l  (p,g   )] + eps ;
pSummary('Startup      Cost'             ,'Clustered-UC') = sum[(p,ta(g)   )                 , pStartupCost (g)                 * vStartUp.l (p,g   )] + eps ;
pSummary('ShutDown     Cost'             ,'Clustered-UC') = sum[(p,ta(g)   )                 , pShutdownCost(g)                 * vShutDown.l(p,g   )] + eps ;
pSummary('ReserveUp    Cost'             ,'Clustered-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResUpCost   * v2ndResUP.l(p,g   )] + eps ;
pSummary('ReserveDw    Cost'             ,'Clustered-UC') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResDwCost   * v2ndResDW.l(p,g   )] + eps ;
pSummary('ENS          Cost'             ,'Clustered-UC') = sum[(sc,p,i    ) $[pScenProb(sc)], pENSCost                         * vENS.l     (sc,p,i)] + eps ;
pSummary('Variable     Cost'             ,'Clustered-UC') = sum[(sc,p,ta(g)) $[pScenProb(sc)], pVariableCost(g)*[vProduct1.l(sc,p,g) + pMinProd(g)* vCommit.l(   p,g)]] + eps ;
pSummary('Curtailment  Cost'             ,'Clustered-UC') = sum[(sc,p,i    ) $[pScenProb(sc)], pWindCurtCost   *[pIntermGen (sc,p,i)               -vIG.l    (sc,p,i)]] + eps ;

* saving of results
pCommitt ('Clustered-UC',   ct,pa(p))                  = vCommit.l  (   p,ct)     + eps ;
pProduct ('Clustered-UC',sc,ct,pa(p)) $[pScenProb(sc)] = vProduct.l (sc,p,ct)*1e3 + eps ;
pENS     ('Clustered-UC',sc,pa(p),i ) $[pScenProb(sc)] = vENS.l     (sc,p,i )*1e3 + eps ;
p2ndResUP('Clustered-UC',ct,pa(p)   )                  = v2ndResUP.l(   p,ct)*1e3 + eps ;
p2ndResDW('Clustered-UC',ct,pa(p)   )                  = v2ndResDW.l(   p,ct)*1e3 + eps ;

* data output to xls file
put TMP put 'par=pCommitt  rdim=2 rng=UC!a1'      / 'par=pProduct  rdim=3 rng=Output!a1'   /
put TMP put 'par=pENS      rdim=3 rng=ENS!a1'     / 'par=pSummary  rdim=1 rng=Summary!a1'  /
put TMP put 'par=p2ndResUP rdim=2 rng=SecResUP!a1'/ 'par=p2ndResDW rdim=2 rng=SecResDW!a1' /
putclose
execute_unload   'tmp.gdx' pProduct pCommitt pENS pSummary p2ndResUP p2ndResDW
execute          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
execute          'del        tmp.gdx                             tmp.txt' ;
*$offtext
*-------------------------------------------------------------------------------
*                    CLUSTERED UC WITH EXTRA CONSTRAINTS
*-------------------------------------------------------------------------------

* activating extra constraints
pExtraClustConst=1;

* reset and new definition of active thermal units
ta(g)          = no  ;
ta(g) $[ct(g)] = yes ;

* activiting ramping constraints for both type of thermal units
tr(g)                  = no  ;
tr(g) $[t(g) or ct(g)] = yes ;

* no startups and shutdowns for individual units
vStartUp.fx (pa(p),t) = 0 ;
vShutDown.fx(pa(p),t) = 0 ;

solve RSDUC using RMIP minimizing vTotalVCost ;
solve  SDUC using  MIP minimizing vTotalVCost ;

pSummary('Final Obj. Function           ','Clustered-UC-WithExtraConst') = SDUC.Objval   + eps ;
pSummary('CPU Time  Model generation [s]','Clustered-UC-WithExtraConst') = SDUC.resGen   + eps ;
pSummary('CPU Time  Model solution   [s]','Clustered-UC-WithExtraConst') = SDUC.resUsd   + eps ;
pSummary('Number of nodes used in MIP   ','Clustered-UC-WithExtraConst') = SDUC.nodUsd   + eps ;
pSummary('Number of variables           ','Clustered-UC-WithExtraConst') = SDUC.numVar   + eps ;
pSummary('Number of discrete variables  ','Clustered-UC-WithExtraConst') = SDUC.numDVar  + eps ;
pSummary('Number of equations           ','Clustered-UC-WithExtraConst') = SDUC.numEqu   + eps ;
pSummary('Number of nonzero elements    ','Clustered-UC-WithExtraConst') = SDUC.numNZ    + eps ;
pSummary('Integrality gap               ','Clustered-UC-WithExtraConst') =abs[SDUC.Objval-RSDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Final relative gap            ','Clustered-UC-WithExtraConst') =abs[SDUC.Objest- SDUC.Objval]/abs[SDUC.Objval] + eps ;
pSummary('Network constraints considered','Clustered-UC-WithExtraConst') = pNetworkConst + eps ;
pSummary('Fix          Cost'             ,'Clustered-UC-WithExtraConst') = sum[(p,ta(g)   )                 , pFixCost     (g)                 * vCommit.l  (p,g   )] + eps ;
pSummary('Startup      Cost'             ,'Clustered-UC-WithExtraConst') = sum[(p,ta(g)   )                 , pStartupCost (g)                 * vStartUp.l (p,g   )] + eps ;
pSummary('ShutDown     Cost'             ,'Clustered-UC-WithExtraConst') = sum[(p,ta(g)   )                 , pShutdownCost(g)                 * vShutDown.l(p,g   )] + eps ;
pSummary('ReserveUp    Cost'             ,'Clustered-UC-WithExtraConst') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResUpCost   * v2ndResUP.l(p,g   )] + eps ;
pSummary('ReserveDw    Cost'             ,'Clustered-UC-WithExtraConst') = sum[(p,ta(g)   )                 , pFixCost     (g)*p2ndResDwCost   * v2ndResDW.l(p,g   )] + eps ;
pSummary('ENS          Cost'             ,'Clustered-UC-WithExtraConst') = sum[(sc,p,i    ) $[pScenProb(sc)], pENSCost                         * vENS.l     (sc,p,i)] + eps ;
pSummary('Variable     Cost'             ,'Clustered-UC-WithExtraConst') = sum[(sc,p,ta(g)) $[pScenProb(sc)], pVariableCost(g)*[vProduct1.l(sc,p,g) + pMinProd(g)* vCommit.l(   p,g)]] + eps ;
pSummary('Curtailment  Cost'             ,'Clustered-UC-WithExtraConst') = sum[(sc,p,i    ) $[pScenProb(sc)], pWindCurtCost   *[pIntermGen (sc,p,i)               -vIG.l    (sc,p,i)]] + eps ;

* saving of results
pCommitt ('Clustered-UC-WithExtraConst',    g,pa(p))                  = vCommit.l  (   p,g )                                + eps ;
pProduct ('Clustered-UC-WithExtraConst',sc,ct,pa(p)) $[pScenProb(sc)] = vProduct.l (sc,p,ct)                           *1e3 + eps ;
pProduct ('Clustered-UC-WithExtraConst',sc, t,pa(p)) $[pScenProb(sc)] =[vProduct1.l(sc,p,t)+pMinProd(t)*vCommit.l(p,t)]*1e3 + eps ;

pENS     ('Clustered-UC-WithExtraConst',sc,pa(p),i ) $[pScenProb(sc)] = vENS.l     (sc,p,i )*1e3 + eps ;
p2ndResUP('Clustered-UC-WithExtraConst',g ,pa(p)   )                  = v2ndResUP.l(   p,g )*1e3 + eps ;
p2ndResDW('Clustered-UC-WithExtraConst',g ,pa(p)   )                  = v2ndResDW.l(   p,g )*1e3 + eps ;

* data output to xls file
put TMP put 'par=pCommitt  rdim=2 rng=UC!a1'        / 'par=pProduct  rdim=3 rng=Output!a1'  /
put TMP put 'par=pENS      rdim=3 rng=ENS!a1'       / 'par=pSummary  rdim=1 rng=Summary!a1' /
put TMP put 'par=p2ndResUP rdim=2 rng=SecResUP!a1'  / 'par=p2ndResDW rdim=2 rng=SecResDW!a1'/
putclose
execute_unload   'tmp.gdx' pProduct pCommitt pENS pSummary p2ndResUP p2ndResDW
execute          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
execute          'del        tmp.gdx                             tmp.txt' ;

$onlisting
