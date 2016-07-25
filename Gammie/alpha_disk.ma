
(** general alpha-disk models **)

(* Charles F. Gammie, for PiTP, 2016 *)

(* physical, mathematical constants *)
pc = 3.08 10^(18)
msun = 2 10^(33)
mp = 1.67 10^(-24)
me = 9.11 10^(-28)
yr = 365.25 24 3600
pi = 3.141592653589793
k = 1.38 10^(-16)
s = 5.67 10^(-5)
A = 7.56 10^(-15)
G = 6.67 10^(-8)
c = 2.998 10^(10)

(* length, mass, accretion rate scalings *)
M = m msun   (* mass *)
Kes = 0.348
Ledd = 4 pi G M c/Kes
NominalEfficiency = 0.1
dMedd = Ledd/(NominalEfficiency c^2)
dM = dm dMedd  (* accretion rate *)
Rs = G M/c^2
R  = x Rs    (* radius *)

(* some alternative scalings *)
(* CVs *)
(*
M = m msun
dM = dm 10^(18)  (* g/s *)
R = x 10^8
*)
(* protoplanetary disks *)
(*
M = m msun
dM = dm 10^(-8) msun/yr
R = x AU
*)

(** opacity regimes **)
(* roster of approximate Rosseland mean opacities, 
   from Bell & Lin *)
(* electron-scattering opacity *)
K0 = Kes
(* bb and ff opacity*)
K1 = 1.5 10^(20) r T^(-5/2)
(* H- *)
K2 = 1.0 10^(-36) r^(1/3) T^(10)
(* Molecules *)
K3 = 10^(-8) r^(2/3) T^3
(* evaporation of metal grains *)
K4 = 2 10^(81) r T^(-24)
(* metal grains *)
K5 = 0.1 T^(1/2)
(* Evaporation of ice grains *)
K6 = 2 10^(16) T^(-7)
(* ice grains *)
K7 = 2 10^(-4) T^2

(* these relations are true in any opacity regime *)
W = (G M/R^3)^(1/2)   
H = cs/W
n = a cs H  (* general version *)
(* this commented version assumed gas pressure dominates *)
(* n = a k T/(mp W) *)
S = 2 H r  (* note: r is density *)
t = S K/2  (* t = midplane optical depth *)
F = (9/8) S n W^2  (* emergent flux, in thermal equilibrium *)
Te = (F/s)^(1/4)  (* surface effective temperature *)

eq1 = T^4 / (3 t Te^4/8)   (* in thermal equilibrium eq1 = 1 *)
eq2 = n S / (dM/(3 pi))    (* in steady state, eq2 = 1 *)

(* zone 0: electron scattering, RP dom *)
Print["zone 0"]
cs := Sqrt[(4/9) A T^4/r]   (* sound speed *)
mu := 0.6 mp                (* mean molecular weight *)
K := K0
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T0 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T0]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r0 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]
T0 = T0 /. r -> r0
S0 = PowerExpand[S /. {r -> r0, T -> T0}]
cs0 = PowerExpand[cs /. {r -> r0, T -> T0}]
H0 = PowerExpand[H /. {r -> r0, T -> T0}]
F0 = PowerExpand[F /. {r -> r0, T -> T0}]
Q0 = cs0 W/(pi G S0)   (* Toomre's Q *)
betar0 = 3 r0 k T0/(mu A T0^4)  (* ratio of gas to radiation pressure *)

	(* outer boundary of zone 0 *)
r01a = PowerExpand[x /. Solve[betar0 == 1,x][[-1]]]
r01b = PowerExpand[x /. Solve[(K0 == K1) /. {r -> r0, T -> T0},x][[-1]]]

(* zone 1: bf,ff opacity, RP dom *)
Print["zone 1"]
cs := Sqrt[(4/9) A T^4/r]
mu := 0.6 mp
K := K1
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T1 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T1]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r1 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]
T1 = T1 /. r -> r1
S1 = PowerExpand[S /. {r -> r1, T -> T1}]
cs1 = PowerExpand[cs /. {r -> r1, T -> T1}]
H1 = PowerExpand[H /. {r -> r1, T -> T1}]
F1 = PowerExpand[F /. {r -> r1, T -> T1}]
tau1 = PowerExpand[ S1 K1 /. {r -> r1, T -> T1}]
Q1 = cs1 W/(pi G S1)
betar1 = 3 r1 k T1/(mu A T1^4)

	(* outer boundary *)
r12a = PowerExpand[x /. Solve[betar1 == 1,x][[-1]]]
r12b = PowerExpand[x /. Solve[PowerExpand[K1/K2 /. 
	{r -> r1, T -> T1}] == 1,x][[-1]]]

(* zone 2: bound-free, free-free, gas pressure dominated *)
Print["zone 2"]
cs := Sqrt[k T/mu]
mu := 0.6 mp
K := K1
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T2 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T2]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r2 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T2 = T2 /. r -> r2
S2 = PowerExpand[S /. {r -> r2, T -> T2}]
cs2 = PowerExpand[cs /. {r -> r2, T -> T2}]
H2 = PowerExpand[H /. {r -> r2, T -> T2}]
F2 = PowerExpand[F /. {r -> r2, T -> T2}]
tau2 = PowerExpand[ S2 K1 /. {r -> r2, T -> T2}]
Q2 = cs2 W/(pi G S2)
betar2 = 3 r2 k T2/(mu A T2^4)

	(* outer boundary *)
r23 = PowerExpand[x /. Solve[PowerExpand[K1/K2 /. 
	{r -> r2, T -> T2}] == 1,x][[-1]]]

(* zone 3: H- opacity, gas pressure dominated *)
Print["zone 3"]
cs := Sqrt[k T/mu]
mu := mp
K := K2
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T3 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T3]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r3 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T3 = T3 /. r -> r3
S3 = PowerExpand[S /. {r -> r3, T -> T3}]
cs3 = PowerExpand[cs /. {r -> r3, T -> T3}]
H3 = PowerExpand[H /. {r -> r3, T -> T3}]
F3 = PowerExpand[F /. {r -> r3, T -> T3}]
tau3 = PowerExpand[ S3 K2 /. {r -> r3, T -> T3}]
Q3 = cs3 W/(pi G S3)
betar3 = 3 r3 k T3/(mu A T3^4)

	(* outer boundary *)
tmp = PowerExpand[K2/K3 /.  {r -> r3, T -> T3}]
r34 = PowerExpand[x /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,x])] == 1, x][[-1]]]

(* zone 4: molecular opacity, gas pressure dominated *)
Print["zone 4"]
cs := Sqrt[k T/mu]
mu := 2.4*mp
K := K3
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T4 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T4]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r4 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T4 = T4 /. r -> r4
S4 = PowerExpand[S /. {r -> r4, T -> T4}]
cs4 = PowerExpand[cs /. {r -> r4, T -> T4}]
H4 = PowerExpand[H /. {r -> r4, T -> T4}]
F4 = PowerExpand[F /. {r -> r4, T -> T4}]
tau4 = PowerExpand[ S4 K3 /. {r -> r4, T -> T4}]
Q4 = cs4 W/(pi G S4)
betar4 = 3 r4 k T4/(mu A T4^4)

	(* outer boundary *)
tmp = PowerExpand[K3/K4 /.  {r -> r4, T -> T4}]
r45 = PowerExpand[x /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,x])] == 1, x][[-1]]]

(* zone 5: metal grain evaporation, gas pressure dom *)
Print["zone 5"]
cs := Sqrt[k T/mu]
mu := 2.4*mp
K := K4
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T5 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T5]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r5 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T5 = T5 /. r -> r5
S5 = PowerExpand[S /. {r -> r5, T -> T5}]
cs5 = PowerExpand[cs /. {r -> r5, T -> T5}]
H5 = PowerExpand[H /. {r -> r5, T -> T5}]
F5 = PowerExpand[F /. {r -> r5, T -> T5}]
tau5 = PowerExpand[ S5 K4 /. {r -> r5, T -> T5}]
Q5 = cs5 W/(pi G S5)
betar5 = 3 r5 k T5/(mu A T5^4)

	(* outer boundary *)
tmp = PowerExpand[K4/K5 /.  {r -> r5, T -> T5}]
r56 = PowerExpand[x /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,x])] == 1, x][[-1]]]


(* zone 6: metal grain opacity, gas pressure dom *)
Print["zone 6"]
cs := Sqrt[k T/mu]
mu := 2.4*mp
K := K5
	(* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T6 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T6]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r6 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T6 = T6 /. r -> r6
S6 = PowerExpand[S /. {r -> r6, T -> T6}]
cs6 = PowerExpand[cs /. {r -> r6, T -> T6}]
H6 = PowerExpand[H /. {r -> r6, T -> T6}]
F6 = PowerExpand[F /. {r -> r6, T -> T6}]
tau6 = PowerExpand[ S6 K5 /. {r -> r6, T -> T6}]
Q6 = cs6 W/(pi G S6)
betar6 = 3 r6 k T6/(mu A T6^4)

	(* outer boundary *)
tmp = PowerExpand[K5/K6 /.  {r -> r6, T -> T6}]
r67 = PowerExpand[x /. Solve[
	PowerExpand[tmp^(1/Exponent[tmp,x])] == 1, x][[-1]]]

(* zone 7: evaporation of ice grains, gas pressure dom *)
Print["zone 7"]
cs := Sqrt[k T/mu]
mu := 2.4*mp
K := K6
        (* solve *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T7 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T7]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r7 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T7 = T7 /. r -> r7
S7 = PowerExpand[S /. {r -> r7, T -> T7}]
cs7 = PowerExpand[cs /. {r -> r7, T -> T7}]
H7 = PowerExpand[H /. {r -> r7, T -> T7}]
F7 = PowerExpand[F /. {r -> r7, T -> T7}]
tau7 = PowerExpand[ S7 K7 /. {r -> r7, T -> T7}]
Q7 = cs7 W/(pi G S7)
betar7 = 3 r7 k T7/(mu A T7^4)

        (* outer boundary *)
tmp = PowerExpand[K6/K7 /.  {r -> r7, T -> T7}]
r78 = PowerExpand[x /. Solve[
        PowerExpand[tmp^(1/Exponent[tmp,x])] == 1, x][[-1]]]

(* zone 8: ice grain opacity, gas pressure dom *)
Print["zone 8"]
cs := Sqrt[k T/mu]
mu := 2.4*mp
K := K7
        (* solve; this one needs to be treated differently
           because of exponents *)
eq1p = PowerExpand[eq1^(1/Exponent[eq1,T])]
T8 = PowerExpand[T /. Solve[eq1p == 1,T][[-1]]]
eq2p = PowerExpand[eq2 /. T -> T8]
eq2p = PowerExpand[eq2p^(1/Exponent[eq2p,r])]
r8 = PowerExpand[r /. Solve[eq2p == 1,r][[-1]]]

T8 = T8 /. r -> r8
S8 = PowerExpand[S /. {r -> r8, T -> T8}]
cs8 = PowerExpand[cs /. {r -> r8, T -> T8}]
H8 = PowerExpand[H /. {r -> r8, T -> T8}]
F8 = PowerExpand[F /. {r -> r8, T -> T8}]
tau8 = PowerExpand[ S8 K8 /. {r -> r8, T -> T8}]
Q8 = cs8 W/(pi G S8)
betar8 = 3 r8 k T8/(mu A T8^4)
        (* outer boundary: boundless *)


(* find Fnu, assuming infinite optically thick, thermal-spectrum disk *)
(* independent of opacity *)
nuV = c/(5500 10^(-8))  (* freq. at center of V band *)
nu = nuV

h = 6.67 10^(-27)       (* Planck *)
F0 = (3/(8 Pi)) W^2 dM  (* integrated surface brightness *)
Teff = (F0/s)^(1/4)
Bnu = 2 h nu^3/c^2/(Exp[h nu/(k Teff)] - 1)

	(* d = distance *)
Fnu = Simplify[ (Cos[i]/d^2) Integrate[2 Pi R Bnu Rs, {x, 0, Infinity}],
	{dm > 0, m > 0} ]

Jy = 10^(-23)
FnuV = 3640 10.^(-0.4 V) Jy  (* V band -> Jy conversion *)

Vexp = PowerExpand[ -48.597 - 1.08574 Log[Fnu] ]


