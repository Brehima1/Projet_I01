# Modélisation et simulation des écoulements de fluides dans la géosphère  
##SUJET DU PROJET : Équation de diffusion non linéaire et solution auto-semblable

---

##  Description
Ce projet traite de la **modélisation et simulation numérique de l’équation de diffusion non linéaire**, avec une application à la **géoscience** (écoulements dans milieux poreux, transferts thermiques, dynamique des plasmas).  
L’objectif est de comparer les comportements **linéaire (m = 1)** et **non linéaire (m > 1)** à travers la **méthode des volumes finis**, et de confronter les résultats numériques à la solution analytique **auto-semblable de Barenblatt**.

---

##  Objectifs
- Étudier l’équation de diffusion non linéaire :  
  \[
  ∂_t u - ∂_x (D(u) ∂_x u) = f(x), \quad D(u) = u^{m-1}, \; m ≥ 1
  \]
- Cas particuliers :  
  - **m = 1** : équation de la chaleur (diffusion linéaire).  
  - **m > 1** : diffusion non linéaire avec support compact (équation des milieux poreux).  
- Implémenter un schéma de **volumes finis implicite**.  
- Utiliser la **méthode de Picard** pour résoudre les non-linéarités.  
- Comparer les solutions numériques avec la solution analytique de **Barenblatt**.  
- Explorer les cas sur :  
  - Domaine infini ℝ.  
  - Domaine semi-infini ℝ⁺.  

---

## Méthodologie
1. **Cas linéaire (m = 1)**  
   - D(u) = 1 → équation de la chaleur.  
   - Validation numérique avec condition initiale de type Dirac.  
   - Comparaison avec la solution analytique gaussienne.  

2. **Cas non linéaire (m > 1)**  
   - Diffusion plus lente avec support compact.  
   - Utilisation d’un schéma implicite (θ = 1).  
   - Itération de Picard pour linéariser le problème à chaque pas de temps.  

3. **Solution de Barenblatt**  
   - Solution auto-semblable de type :  
     \[
     u(x,t) \sim t^{-α} F\left(\frac{x}{t^β}\right)
     \]  
   - Comparaison des résultats numériques et analytiques pour plusieurs valeurs de m.  

4. **Implémentation**  
   - `darcy_1D.py` : construction de la matrice et second membre.  
   - `heat_1D.py` : méthode de Picard pour la diffusion non linéaire.  
   - Notebooks Python : expériences numériques sur ℝ et ℝ⁺.  
   - `barenblatt_solution.py` : solution exacte de Barenblatt pour comparaison.  

---

##  Résultats
- **m = 1 (cas linéaire)** :  
  - La solution numérique coïncide avec la gaussienne attendue.  
  - Validation du schéma de volumes finis.  

- **m > 1 (cas non linéaire)** :  
  - La solution a un **support compact** → diffusion plus localisée.  
  - Pour m grand, la propagation ralentit fortement, créant des fronts marqués.  
  - Les solutions numériques suivent fidèlement la solution de Barenblatt.  

- **Comparaison linéaire vs non linéaire** :  
  - m → 1 : la solution de Barenblatt tend vers la gaussienne.  
  - m > 1 : apparition d’un profil à support borné, contraste avec le cas linéaire.  

---

##  Compétences démontrées
- Analyse et résolution d’**équations aux dérivées partielles non linéaires**.  
- Implémentation d’un **schéma de volumes finis implicite**.  
- Utilisation de la **méthode de Picard** pour gérer les non-linéarités.  
- Validation numérique par confrontation avec une solution analytique (Barenblatt).  
- Programmation scientifique en **Python (NumPy, SciPy, Matplotlib, Jupyter)**.  

---

##  Organisation
- `Projet_I01_2024_2025.pdf` → Rapport complet avec démonstrations et résultats.  
- `code/` → Scripts Python (`darcy_1D.py`, `heat_1D.py`, `barenblatt_solution.py`).  
- `notebooks/` → Expériences numériques (diffusion linéaire et non linéaire).  
- `figures/` → Comparaisons des solutions analytiques et numériques.  
- `README.md` → Présentation claire et professionnelle du projet.  

---

##  Perspectives
- Utilisation de schémas numériques **d’ordre supérieur** et maillages adaptatifs.  
- Extension à des **géométries plus complexes** et **dimensions supérieures**.  
- Application directe aux **écoulements souterrains, milieux poreux, thermique**.  
- Intégration avec des solveurs avancés (FEniCS, Firedrake).  

---

##  Auteur
**Bréhima Samaké**  
2025  
 Projet académique – Modélisation et simulation des écoulements de fluides  

---

## Conlusion : 
- Maîtrise de la **modélisation physique et mathématique de phénomènes complexes**.  
- Compétences avancées en **analyse numérique, méthodes implicites et résolution non linéaire**.  
- Capacité à **implémenter et valider des solutions numériques** rigoureuses.  
- Expérience directement transférable en **R&D, ingénierie, data science et simulation physique**.  
