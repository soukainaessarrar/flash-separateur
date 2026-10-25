"""
================================================================================
  FLASH ISOTHERME — SEPARATEUR DIPHASIQUE
  Equations MESH + Rachford-Rice (loi de Raoult)
================================================================================

Schema :
                    F, zA, zB, T, P
                         |
                         v
                   +-----------+
                   | SEPARATEUR|   (flash isotherme)
                   |   MESH    |
                   +-----+-----+
                         |
            +------------+------------+
            |                         |
            v                         v
     VAPEUR (V, yA, yB)        LIQUIDE (L, xA, xB)

Entrees  : xA, xB (fractions molaire alimentation), T (degC), P (mmHg)
Sorties  : compositions vapeur et liquide + debits

Melange  : Benzene (A) / Toluene (B)
Antoine  : log10(Psat[mmHg]) = A - B/(C + T[degC])
Equilibre: Ki = Psat_i(T) / P
Bilan    : F*zi = V*yi + L*xi   avec  beta = V/F
           yi = Ki * xi
           sum xi = 1 , sum yi = 1
           -> equation de Rachford-Rice sur beta

Lancer :
  python flash_separateur.py

Exemple (T = 80 degC) :
  xA = 0.5
  xB = 0.5
  T  = 80
  P  = 500
================================================================================
"""

from pathlib import Path
import numpy as np

# ----------------------------------------------------------------------
# Donnees composants
# ----------------------------------------------------------------------
NOM_A, NOM_B = "Benzene", "Toluene"
MA, MB = 78.11, 92.14  # g/mol
ANTOINE_A = (6.90565, 1211.033, 220.790)
ANTOINE_B = (6.95464, 1344.800, 219.482)
F = 100.0  # debit alimentation (kmol/h)

OUT_DIR = Path(__file__).resolve().parent


def psat(coef, T):
    A, B, C = coef
    return 10 ** (A - B / (C + T))


def K_values(T, P):
    """Equation E (equilibre) : Ki = Psat_i / P"""
    return np.array([psat(ANTOINE_A, T), psat(ANTOINE_B, T)]) / P


def f_rr(beta, z, K):
    """Fonction de Rachford-Rice f(beta) = 0"""
    return float(np.sum(z * (K - 1.0) / (1.0 + beta * (K - 1.0))))


def resoudre_beta(z, K, tol=1e-12):
    """Resout f(beta)=0 (bissection). Gere liquide seul / vapeur seule."""
    if f_rr(0.0, z, K) <= 0:
        return 0.0, "1 phase : LIQUIDE seulement (sous le point de bulle)"
    if f_rr(1.0, z, K) >= 0:
        return 1.0, "1 phase : VAPEUR seulement (au-dessus du point de rosee)"
    a, b = 0.0, 1.0
    while b - a > tol:
        m = 0.5 * (a + b)
        if f_rr(a, z, K) * f_rr(m, z, K) <= 0:
            b = m
        else:
            a = m
    return 0.5 * (a + b), "2 phases : VAPEUR + LIQUIDE (separateur diphasique)"


def flash_separateur(zA, zB, T, P):
    """
    Flash TP du separateur (MESH sans H : T et P imposes).

    M : F*zi = V*yi + L*xi
    E : yi = Ki*xi
    S : sum x = 1, sum y = 1
    """
    s = zA + zB
    if s <= 0:
        raise ValueError("xA + xB doit etre > 0")
    zA, zB = zA / s, zB / s
    z = np.array([zA, zB])
    K = K_values(T, P)
    beta, etat = resoudre_beta(z, K)
    x = z / (1.0 + beta * (K - 1.0))  # liquide sortant
    y = K * x                         # vapeur sortante
    return z, K, beta, x, y, etat


def lire_float(invite, defaut=None):
    while True:
        brut = input(invite).replace(",", ".").strip()
        if brut == "" and defaut is not None:
            return float(defaut)
        try:
            return float(brut)
        except ValueError:
            print("  -> nombre invalide. Exemple : 0.5 ou 80")


def afficher_et_sauver(z, K, beta, x, y, etat, T, P):
    V, L = beta * F, (1.0 - beta) * F
    M = np.array([MA, MB])
    wx = x * M / np.sum(x * M)
    wy = y * M / np.sum(y * M)
    bilan = F * z - (V * y + L * x)

    lignes = []
    def out(s=""):
        print(s)
        lignes.append(s)

    out("=" * 66)
    out(f" SEPARATEUR FLASH  {NOM_A} (A) / {NOM_B} (B)")
    out(" Methode : MESH + Rachford-Rice (Raoult)")
    out("=" * 66)
    out("")
    out("  Alimentation F -----> [ SEPARATEUR ] -----> Vapeur  V (yA, yB)")
    out("                        [   T, P     ] -----> Liquide L (xA, xB)")
    out("")
    out(f" ENTREES")
    out(f"   zA = {z[0]:.4f} | zB = {z[1]:.4f} | F = {F} kmol/h")
    out(f"   T  = {T} degC  | P  = {P} mmHg")
    out(f"   Psat A = {psat(ANTOINE_A, T):7.1f} mmHg | Psat B = {psat(ANTOINE_B, T):7.1f} mmHg")
    out(f"   K_A = {K[0]:.4f} | K_B = {K[1]:.4f}")
    out(f"   Etat du separateur : {etat}")
    out("-" * 66)
    out(f" DEBITS (sortie du separateur)")
    out(f"   beta = V/F = {beta:.5f}")
    out(f"   Vapeur  V = {V:.2f} kmol/h")
    out(f"   Liquide L = {L:.2f} kmol/h")
    out("-" * 66)
    out(" COMPOSITIONS — PHASE VAPEUR")
    out(f"   yA = {y[0]:.4f} | yB = {y[1]:.4f} | somme = {y.sum():.4f}")
    out(f"   wA = {wy[0]:.4f} | wB = {wy[1]:.4f}  (massique)")
    out("-" * 66)
    out(" COMPOSITIONS — PHASE LIQUIDE")
    out(f"   xA = {x[0]:.4f} | xB = {x[1]:.4f} | somme = {x.sum():.4f}")
    out(f"   wA = {wx[0]:.4f} | wB = {wx[1]:.4f}  (massique)")
    out("-" * 66)
    out(f" Verification bilan M : F*z - (V*y + L*x) = {np.round(bilan, 12)}")
    out("=" * 66)

    rapport = OUT_DIR / "resultat_separateur.txt"
    rapport.write_text("\n".join(lignes) + "\n", encoding="utf-8")
    print(f"\n Resultats enregistres : {rapport}")


def main():
    print("=" * 66)
    print(f" FLASH / SEPARATEUR  {NOM_A} (A) / {NOM_B} (B)")
    print(" Entrez xA, xB, T, P  (T typique = 80 degC)")
    print("=" * 66)
    print()

    xA = lire_float("  xA (alimentation) = ")
    xB = lire_float("  xB (alimentation) = ")
    T = lire_float("  T  (degC) [defaut 80] = ", defaut=80)
    P = lire_float("  P  (mmHg)           = ")

    if abs(xA + xB - 1.0) > 1e-6:
        print(f"\n  Note : xA+xB = {xA + xB:.4f} != 1 -> normalisation.")

    z, K, beta, x, y, etat = flash_separateur(xA, xB, T, P)
    print()
    afficher_et_sauver(z, K, beta, x, y, etat, T, P)


if __name__ == "__main__":
    main()
