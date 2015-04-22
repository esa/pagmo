/*****************************************************************************
 *   Copyright (C) 2015 The PaGMO development team,                          *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef PAGMO_ALGORITHM_WORHP_SERIALIZATION_H
#define PAGMO_ALGORITHM_WORHP_SERIALIZATION_H

#include <worhp/worhp.h>

// This is from WORHP 1.8.0

template <class Archive>
void serialize(Archive &ar, QPParamsStruct &cb, const unsigned int)
{
  ar & cb.ipBarrier;
  ar & cb.ipComTol;
  ar & cb.ipFracBound;
  ar & cb.ipMinAlpha;
  ar & cb.ipRelaxDiv;
  ar & cb.ipRelaxMax;
  ar & cb.ipRelaxMin;
  ar & cb.ipRelaxMult;
  ar & cb.ipResTol;
  ar & cb.lsTol;
  ar & cb.nsnBeta;
  ar & cb.nsnKKT;
  ar & cb.nsnMinAlpha;
  ar & cb.nsnSigma;

  ar & cb.ipLsMethod;
  ar & cb.lsItMaxIter;
  ar & cb.lsItMethod;
  ar & cb.lsItPrecondMethod;
  ar & cb.lsRefineMaxIter;
  ar & cb.maxIter;
  ar & cb.method;
  ar & cb.nsnLsMethod;
  ar & cb.printLevel;

  ar & cb.ipTryRelax;
  ar & cb.lsScale;
  ar & cb.lsTrySimple;
  ar & cb.nsnGradStep;
  ar & cb.scaleIntern;
  ar & cb.strict;
}

template <class Archive>
void serialize(Archive &ar, ParamsStruct &cb, const unsigned int) 
{
/* Pointer to QP parameter structure */
  ar & cb.qp;

 /* Debug marker. Used to find memory alignment/padding issues */
  ar & cb.DebugMarker05;

 /*------------------------------------------------
  *  Real parameters
  *------------------------------------------------ */

 /* Tolerance for acceptable feasibility */
  ar & cb.AcceptTolFeas;

 /* Tolerance for acceptable optimality */
  ar & cb.AcceptTolOpti;

 /* Trial stepsize decrease factor for Armijo rule */
  ar & cb.ArmijoBeta;

 /* Initial alpha for Armijo rule */
  ar & cb.ArmijoMaxAlpha;

 /* Lower bound on alpha for Armijo rule */
  ar & cb.ArmijoMinAlpha;

 /* Lower bound on alpha for Armijo rule during recovery */
  ar & cb.ArmijoMinAlphaRec;

 /* Scale factor for linearised descent check in Armijo rule */
  ar & cb.ArmijoSigma;

 /* Update factor for Betts' Hessian regularisation */
  ar & cb.BettsFactor;

 /* Smallest eigenvalue of the regularised Hessian */
  ar & cb.BettsPoint;

 /* Factor in determining active constraints by KKT */
  ar & cb.BoundTolFac;

 /* (experimental) */
  ar & cb.CorStepBettsSum;

 /* (experimental) */
  ar & cb.CorStepConvio;

 /* (experimental) */
  ar & cb.CorStepConStop;

 /* (experimental) */
  ar & cb.CorStepPFactor;

 /* (experimental) */
  ar & cb.CorStepPMax;

 /* Upper bound used by Fritz-John heuristic */
  ar & cb.CheckFJ;

 /* Block BFGS curvature condition bound */
  ar & cb.CurvBCond;

 /* Block BFGS curvature condition regularisation factor */
  ar & cb.CurvBFac;

 /* BFGS Curvature condition bound */
  ar & cb.CurvCond;

 /* BFGS curvature condition regularisation factor */
  ar & cb.CurvFac;

 /* Feasibility tolerance for no-objective feasible mode */
  ar & cb.FeasibleInitTol;

 /* Finite difference perturbation */
  ar & cb.FidifEps;

 /* Factor in Focus-on-Feasibility mode */
  ar & cb.FocusOnFeasFactor;

 /* Upper bound for numbers to be regarded as finite */
  ar & cb.Infty;

 /* Tolerance for unboundedness detection heuristic */
  ar & cb.InftyUnbounded;

 /* IP complementarity tolerance in initial multiplier estimate */
  ar & cb.LMestQPipComTol;

 /* IP residual tolerance in initial multiplier estimate */
  ar & cb.LMestQPipResTol;

 /* Lowpass-filter update factor for objective values */
  ar & cb.LowPassAlphaF;

 /* Lowpass-filter update factor for constraint values */
  ar & cb.LowPassAlphaG;

 /* Lowpass-filter update factor for merit function values */
  ar & cb.LowPassAlphaMerit;

 /* Any pivot whose modulus is less than this is treated as zero by MA97 */
  ar & cb.MA97small;

 /* Relative pivot tolerance of MA97 */
  ar & cb.MA97u;

 /* Threshold of meritfunction gradient for increasing Hessian regularisation */
  ar & cb.MeritGradTol;

 /* Penalty update parameter factor for MeritFunction = 4 */
  ar & cb.PenUpdEpsKFac;

 /* Penalty update parameter factor for MeritFunction = 3 */
  ar & cb.PenUpdEpsBar;

 /* Max penalty for MeritFunction = 4 */
  ar & cb.PenUpdMaxDeltaK;

 /* Max factor for increasing penalty for MeritFunction = 4 */
  ar & cb.PenUpdMaxFac;

 /* Penalty update parameter for MeritFunction = 3 */
  ar & cb.PenUpdRBar;

 /* (currently unused) Relative precision of objective */
  ar & cb.PrecisionF;

 /* (currently unused) Relative precision of constraints */
  ar & cb.PrecisionG;

 /* (currently unused) Scaling factor for QP */
  ar & cb.QPscaleParam;

 /* Upper bound for accepting the constraint relaxation variable */
  ar & cb.RelaxMaxDelta;

 /* Upper bound on the constraint relaxation penalty */
  ar & cb.RelaxMaxPen;

 /* Update factor for the constraint relaxation penalty */
  ar & cb.RelaxRho;

 /* Initial value of the constraint relaxation penalty */
  ar & cb.RelaxStart;

 /* Value to scale large objective functions to */
  ar & cb.ScaleFacObj;

 /* Upper bound on resulting matrix norm for QP scaling */
  ar & cb.ScaleFacQP;

 /* Initial value for Betts' update dampening term */
  ar & cb.StartBettsTau;

 /* Timeout in seconds */
  ar & cb.Timeout;

 /* Complementarity tolerance */
  ar & cb.TolComp;

 /* Feasibility tolerance */
  ar & cb.TolFeas;

 /* Optimality tolerance */
  ar & cb.TolOpti;

 /* (experimental) */
  ar & cb.TolWeakActive;

 /* Upper bound on constraint violation for too-big heuristic */
  ar & cb.TooBigCV;

 /* Upper bound on KKT values for too-big heuristic */
  ar & cb.TooBigKKT;

 /* Machine epsilon */
  ar & cb.eps;

 /* Increase factor for estimated integer workspace requirement */
  ar & cb.IncreaseIWS;

 /* Increase factor for estimated real workspace requirement */
  ar & cb.IncreaseRWS;

 /* Constraint violation decrease factor in Filter acceptance check */
  ar & cb.FilterGammaCV;

 /* Objective decrease factor in Filter acceptance check */
  ar & cb.FilterGammaF;

 /* Safety factor for alphamin calculation by Filter */
  ar & cb.GammaAlpha;

 /* Increase factor for Betts' update dampening term */
  ar & cb.IncBettsTau;

 /* Larger increase factor for Betts' update dampening term */
  ar & cb.IncBettsTauMore;

 /* Lower bound for Betts' update dampening term */
  ar & cb.MinBettsTau;

 /* Decrease factor for Betts' update dampening term */
  ar & cb.ReduceBettsTau;

 /* Start tolerance for successful termination of iterative refinement due to perturbation in constraints */
  ar & cb.RefineStartTol;

 /* Maximum allowed relaxation to apply feasibility refinement */
  ar & cb.RefineMaxRelax;

 /* Maximum allowed regularisation of the hessian CAUTION absolute value */
  ar & cb.RefineMaxHMReg;

 /* Filter switching condition parameter */
  ar & cb.SwitchingDelta;

 /* Filter switching condition parameter */
  ar & cb.SwitchingSF;

 /* Filter switching condition parameter */
  ar & cb.SwitchingSCV;

 /*------------------------------------------------
  *  Integer parameters
  *------------------------------------------------ */

 /* Armijo recovery strategies */
  ar & cb.Ares[NAres];

 /* Choose BFGS method (0: dense, 1-3: block, 100+: sparse) */
  ar & cb.BFGSmethod;

 /* Restart BFGS update after this many iterations */
  ar & cb.BFGSrestart;

 /* Block size parameter used by certain BFGS methods */
  ar & cb.BFGSmaxblockSize;

 /* Block size parameter used by certain BFGS methods */
  ar & cb.BFGSminblockSize;

 /* (experimental) */
  ar & cb.CorStepMaxIter;

 /* (experimental) */
  ar & cb.CorStepMethod;

 /* (experimental) */
  ar & cb.CorStepMode;

 /* Select method to determine graph colouring groups */
  ar & cb.GroupMethod;

 /* Enable XML logfiles and writing interval */
  ar & cb.LogLevel;

 /* Enable XML result logging and detail level */
  ar & cb.LogResult;

 /* Ordering used by MA97 */
  ar & cb.MA97ordering;

 /* Scaling used by MA97 */
  ar & cb.MA97scaling;

 /* Print level used by MA97 */
  ar & cb.MA97print;

 /* Node amalgation, controls merging in elimination tree by MA97 */
  ar & cb.MA97nemin;

 /* Upper bound to Reverse Communication calls */
  ar & cb.MaxCalls;

 /* Maximum number of Force recovery strategy steps */
  ar & cb.MaxForce;

 /* (experimental) */
  ar & cb.MaxGPart;

 /* Upper bound on major iterations */
  ar & cb.MaxIter;

 /* Select merit function and penalty update [0, 3..5] */
  ar & cb.MeritFunction;

 /* Select (1) Meritfunction or (3) Filter globalisation */
  ar & cb.NLPmethod;

 /* NLP print level [-1..4] */
  ar & cb.NLPprint;

 /* Select method to determine graph colouring pairgroups */
  ar & cb.PairMethod;

 /* Penalty update parameter */
  ar & cb.PenUpdEpsKSequence;

 /* 0 - Deactivated, 1 - After first feasible iterate, 2 - Always on, Activates iterative refinement due to perturbation in constraints using parametric sensitivities */
  ar & cb.RefineFeasibility;

 /* Enable automatic Hessian structure generation or checking */
  ar & cb.UserHMstructure;

 /* Control activation of Filter acceleration heuristics */
  ar & cb.MaxLScounter;

 /* Select Hessian regularisation strategy in Filter */
  ar & cb.RegStrategy;

 /*------------------------------------------------
  *  Logical parameters
  *------------------------------------------------ */

 /* Enable automatic QP recovery */
  ar & cb.AutoQPRecovery;

 /* Enable structural checking of DF */
  ar & cb.CheckStructureDF;

 /* Enable structural checking of DG */
  ar & cb.CheckStructureDG;

 /* Enable structural checking of HM */
  ar & cb.CheckStructureHM;

 /* (experimental) */
  ar & cb.CorStepRecoveryDX;

 /* F and G cannot be evaluated separately */
  ar & cb.FGtogether;

 /* Enable Fritz-John and non-differentiable check heuristics */
  ar & cb.FJandND;

 /* Activate dual feasibility mode */
  ar & cb.FeasibleDual;

 /* Activate initial feasibility mode */
  ar & cb.FeasibleInit;

 /* Activate feasible-only mode */
  ar & cb.FeasibleOnly;

 /* Approximate Hessian by finite differences (otherwise BFGS) */
  ar & cb.FidifHM;

 /* Use central finite difference quotient for first derivatives */
  ar & cb.FirstDifCentral;

 /* Enable Focus-on-Feasibility mode */
  ar & cb.FocusOnFeas;

 /* Enable initial Lagrange multiplier estimate */
  ar & cb.InitialLMest;

 /* Counter for changed parameters. Internal use only. */
  ar & cb.internalParChanged;

 /* Save acceptable solutions as fallback */
  ar & cb.KeepAcceptableSol;

 /* Control Lagrange multiplier update */
  ar & cb.LinMult;

 /* Enable lowpass-filter termination criterion */
  ar & cb.LowPassFilter;

 /* Use BLAS level 3 (dgemm) in MA97 */
  ar & cb.MA97blas3;

 /* Use multifrontal-style forward solve of MA97 */
  ar & cb.MA97mf;

 /* Not to be included into a parameter file! */
  ar & cb.MatrixCC;

 /* Introduce one relaxation variable for every constraint */
  ar & cb.MoreRelax;

 /* Not to be included into a parameter file! */
  ar & cb.QuadraticProblem;

 /* Use Steffensen Extrapolation during Feasibility Refinement */
  ar & cb.SteffensenOnRefine;

 /* Activates new iterative refinement of constraints only when Armijo alpha equals one */
  ar & cb.RefineOnlyOnAlpha;

 /* Do restoration until a feasible solution is found */
  ar & cb.RestUntilFeas;

 /* Scale constraints in every iteration */
  ar & cb.ScaleConIter;

 /* Use a scaled perturbation for finite differences */
  ar & cb.ScaledFD;

 /* Scale KKT conditions */
  ar & cb.ScaledKKT;

 /* Scale the objective function */
  ar & cb.ScaledObj;

 /* Scale some matrices handed to the QP */
  ar & cb.ScaledQP;

 /* Evaluate QP search direction regardless of convergence */
  ar & cb.TakeQPSol;

 /* Enable too-big termination heuristics */
  ar & cb.TooBig;

 /* Activates update of lagrange multipliers during correction step */
  ar & cb.UpdateMu;

 /* Objective gradient values supplied by caller */
  ar & cb.UserDF;

 /* Jacobian values supplied by caller */
  ar & cb.UserDG;

 /* Hessian values supplied by caller */
  ar & cb.UserHM;

 /* Hessian values supplied by caller */
  ar & cb.UserZenDGp;

 /* Hessian values supplied by caller */
  ar & cb.UserZenDLxp;

 /* Gradient values supplied by caller */
  ar & cb.UserZenDLp;

 /* Hessian values supplied by caller */
  ar & cb.UserZenDLpp;

 /* Run Zen module after successful termination */
  ar & cb.UseZen;

 /* Check maximum of secure perturbation when updating solution */
  ar & cb.ZenCheckMaxPert;

 /* false: use LU from last QP step; true: renew LU decomposition. */
  ar & cb.ZenRenewLU;
  ar & cb.ZenFDnewMethod;

 /* (experimental) */
  ar & cb.WeakActiveSet;

 /* Use a constant lower bound on Armijo stepsize in Filter */
  ar & cb.AlphaMinConst;

 /* Activate accelerating heuristics for Filter */
  ar & cb.IgnoreFilterCrit;

 /* Filter heuristic to save Armijo iterations */
  ar & cb.FilterBisecAlpha;

 /* Filter heuristic to save Armijo iterations */
  ar & cb.FilterIntersecAlpha;

 /* Select max-norm instead of 1-norm in Filter */
  ar & cb.MaxNorm;

 /* Enables Filter-reinitialisation accelerating heuristic */
  ar & cb.ReinitFilter;

 /* Debug marker. Used to find memory alignment/padding issues */
  ar & cb.DebugMarker06;

 /* Automatically added initialisation flag.  */
  ar & cb.initialised;
}

#endif // PAGMO_ALGORITHM_WORHP_SERIALIZATION_H
