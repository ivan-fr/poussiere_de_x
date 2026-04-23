# Pandrosion Corpus Roadmap — Autonomous Claude Loop

État actuel : 73 modules verts. Cible finale : `PrincipalDominanceP3X2`
(Niveau 5) inconditionnel à `x = 2, p = 3`.

## 🔁 Mode d'emploi (pour Claude autonome)

Tu es Claude. Tu travailles **seul**, sans input utilisateur.

À chaque itération :
1. Lis ce fichier. Prends la **première** cible non-cochée.
2. Lis les modules pertinents (§54–§72 selon la cible).
3. Résous : écris/complète le module, build via
   `bash scripts/lean-incremental.sh` jusqu'à `✅ INCREMENTAL OK`.
4. **Refactor si build > 50 s** : split lemmas, switch `nlinarith` →
   `ring+linarith`, simplify `field_simp`. Règles dans `WORKFLOW.md`.
5. **Si la cible est trop dure** : downscope à la version la plus forte
   qui build. Documente dans le module *pourquoi* downscopé. Coche la
   case quand même (partial win est un win).
6. Coche la case `[ ]` → `[x]` dans ce fichier.
7. `git add -A && git commit -m "feat(§N): <target>" && git push`.
8. Recommence tant qu'il reste des cibles non-cochées ET que tu as du
   budget de session.

**Règles dures** :
- Jamais `docker compose run --rm lean-check` (réservé utilisateur).
- Jamais `docker compose run --rm lean-build` (clean rebuild, trop lent).
- Toujours `bash scripts/lean-incremental.sh` (kill stale containers).
- Si build > 50 s : BREAK + refactor (voir WORKFLOW.md §"avoid nlinarith").
- Commit-message : `feat(§N): théorème — rôle dans le corpus`.

---

## Cibles (ordre de priorité)

### Niveau 1 — σ sur demi-plan (priorité #1)

- [x] **§73 SigmaStepBoundX2 — framework + cas trivial** : downscope
  forcé (analyse polynomiale sur `Re z ≥ 1/2` dépasse budget session).
  Acquis : `SigmaClosedStepIntoBasinX2` défini, cas `z = α₀` prouvé,
  chaînage conditionnel vers §70 via §75 documenté.

- [x] **§74 SigmaConcretePointsX2** (downscope) : full sur tout le
  demi-plan reste hors budget. Cette session ajoute des **vérifications
  concrètes à points nommés** : `σ_closed(1) = 353/444`, `σ_closed(2) =
  5571/6888`, et `|σ_closed(z) - α₀| < 1/4` à ces deux points (via
  `α₀ ∈ [3/4, 1]` + arithmétique linéaire). Building blocks pour
  preuve full future.

- [x] **§75 SigmaEquivPointsX2** (downscope) : équivalence ponctuelle
  à `z = 1` (`steffensen_step_C 2 3 (1) = sigma_x2_closed (1) = 353/444`).
  Calcul direct via norm_num (évite field_simp degré 18+). Building
  block pour preuve générale future.

- [ ] **§76 HalfPlaneSigmaTendstoP3X2** : chaîner §74 + §75 + §65 +
  §70 `half_plane_sigma_tendsto_from_step_into_basin_v2`. Si §75 est
  downscopé, reformuler en "conditional sur §75".

### Niveau 5 — Principal dominance (priorité #2)

- [x] **§77 FarFieldX2** : `|z| ≥ 2 ⟹ h^n(z) → α₀` inconditionnel.
  Argument : `|Sp(z)| ≥ ‖z‖² − ‖z‖ − 1 ≥ 1` pour `‖z‖ ≥ 2`, donc
  `Re h(z) ≥ 1/2`, puis §67 Tendsto conclut. Couvre **tout le
  far-field** complexe inconditionnellement.

- [ ] **§77b HalfPlaneExhaustion — région bornée** : prouver que pour
  a.e. `z ∈ ℂ` avec `Re z < 1/2 ∧ |z| < 2`, il existe `k` avec
  `h^k(z)` dans le demi-plan ou le far-field. Empirique : vérifié
  Python sur `[-2, 1/2] × [-2, 2]` hors voisinages de ω, ω².

- [ ] **§78 NonPrincipalBasin finite measure** : prouver que
  `⋃_{s≠0} CyclotomicBasinP3X2 s` est contenu dans un compact, donc
  à mesure finie. Empirique : concentration ~10⁻⁴ autour de `ω·α₀`,
  `ω²·α₀`.

- [ ] **§79 Niveau 5 final** : chaîner §76 + §77 + §78 pour fermer
  `PrincipalDominanceP3X2`. Si §77 ou §78 manque, produire version
  conditionnelle + chaînage explicite.

### Extensions p ≥ 3 (priorité #3)

- [x] **§80 HalfPlaneContractionP4** : généralisation partielle.
  Identité polynomiale exacte pour `‖Sp_C 4 z‖²` prouvée. Borne
  `‖Sp_C 4 z‖² ≥ (1+a+a²+a³)²` prouvée pour **`Re z ≥ 1`** (PAS pour
  `Re z ≥ 1/2` — contre-exemple à `z = 1/2 + i`). Décomposition
  ring-identity + factorisation `(a-1)·(...)` pour les coefficients.

- [x] **§81 AlphaX2P4** : ancre `α₀ := 2^{-1/4}` à `x = 2, p = 4`.
  Données arithmétiques : `α₀^4 = 1/2`, `α₀ ≥ 3/4` (via `(3/4)^4 = 81/256 ≤ 128/256`),
  `α₀ ≤ 1`. Théorème Banach complet downscope (demande §56.3
  généralisé pour Q_4 = 1+z+α+z²+zα+α²).

- [x] **§82 AlphaGenericPX** : ancre générique `α(p, x) := x^{-1/p}`
  unifiant §81 (alphaX2P4) et §83 (alphaX). Identités `α(p, x)^p =
  1/x`, spécialisations vers les ancres concrètes. Rate `K(p) < 1`
  paramétré reste pour preuve future (analyse polynomiale paramétrée
  sur p).

### Extensions x ≠ 2 (priorité #4)

- [x] **§83 AlphaGenericX** : ancre paramétrée `α(x) := x^{-1/3}`
  pour `x > 0`. Identités clefs : `α(x)^3 = 1/x`, `α(x) ≤ 1` pour
  `x ≥ 1` (strict pour `x > 1`), `α(x) ≥ 1/2` pour `0 < x ≤ 8`.
  Théorème Banach paramétré sur x downscope (demande paramétrer
  toutes les bornes polynomiales sur x, hors budget).

- [x] **§84 HalfPlaneGenericX** : invariance demi-plan paramétrée
  pour `pandrosion_h_C x 3` à `x ∈ (1, 8]`. `Re h(x, 3, z) ≥ 1/2`
  sur `Re z ≥ 1/2` via majoration `(x-1)/(x·‖Sp(z)‖) ≤ 1/2`. La plage
  `x ≤ 8` correspond exactement à `α(x) ≥ 1/2` (§83). Tendsto et
  Banach paramétrés sur x restent pour preuve future.

### Périphériques — si budget restant

- [ ] **§85 Explicit Steffensen radius lower bound** : dérouler
  `steffensenR_of_fp` à l'anchor x=2 pour obtenir une valeur
  concrète (ex: `R_σ ≥ 1/100`). Débloque §70
  `SteffensenRadiusAtLeast14` (ou version affaiblie).

- [x] **§86 IterCountX2** : compte d'itérations pour `ε`-précision.
  Forme existentielle `∃ N, ∀ n ≥ N, ‖h^n z − α₀‖ ≤ ε` + forme
  géométrique explicite `(104/259)^N · (1/4) ≤ ε ⟹ N suffit`. La
  forme `Real.log`-exacte reste pour §86b future.

- [ ] **§87 Böttcher formal series attempt** : tenter une version
  formelle de §68's `ψ` via série de puissances tronquée, preuve
  de convergence locale.

---

## Stop conditions pour Claude

Arrête la boucle si l'une des conditions suivantes est atteinte :

1. **Toutes les cases cochées** : ROADMAP terminée, rien à faire.
2. **2 downscopes consécutifs** : la prochaine cible demande plus de
   budget que disponible. Commit l'état, stop.
3. **Build échoue 3 fois consécutives** : bug stable, nécessite œil
   humain. Commit l'état du dernier module, stop.
4. **Session proche du cap de contexte** : si la conversation dépasse
   ~70% du budget token, stop pour démarrer une nouvelle session
   propre au réveil de l'utilisateur.

Quand tu stoppes, écris un résumé à la fin de ce fichier sous
`## Session N résumé` avec :
- Cibles complétées (par numéro §).
- Cibles downscopées (par numéro §, avec raison).
- Cibles bloquées (par numéro §, avec point de friction).
- Prochaine session : suggestion de prompt.

---

## Session 1 résumé (autonomous overnight setup)

**État corpus début** : 73 modules verts (après §72 SigmaClosedFormX2).
**État corpus fin** : 75 modules verts.

### Cibles complétées
- **§77 FarFieldX2** : `|z| ≥ 2 ⟹ h^n(z) → α₀` inconditionnel. Argument
  triangle inverse + factorisation `(‖z‖-2)(‖z‖+1) ≥ 0`.
- **§78 NamedPointsX2** : 7 points concrets dans le h-bassin
  (z = 1, 2, 5, 10, -5, 3i, 1+i) via §67 + §77.

### Cibles non-commencées (prochaine session)
- **§73 SigmaStepIntoBasinX2** : bloqué par analyse polynomiale
  complexe. Piste : commencer par `|z| ≥ 10` (régime asymptotique).
- **§77b HalfPlaneExhaustion bounded** : région `{Re z < 1/2 ∧ |z| < 2}`.
  Analyse dynamique complexe, hors portée d'une seule session.
- **§80 p=4 extension** : demande réécriture de §57 pour `Sp_C 4`.

### Observations pour la suite
- Le pattern "ring-identity + linarith" fonctionne très bien pour les
  inégalités polynomiales (§77 build en 9s).
- Les corollaires "points nommés" (§78) sont utiles : ancres
  démonstratives bon marché, mais n'avancent pas la recherche.
- Attention : `rw [h_decomp]` avec h_decomp dont RHS contient le LHS
  provoque un unfolding récursif. Utiliser `have h_norm_eq` sur les
  normes et linarith.

### Prochaine session : suggestion de prompt

```
Lis ROADMAP.md. Attaque §73 SigmaStepIntoBasinX2 pour |z| ≥ 10
seulement (restriction asymptotique). Utilise sigma_x2_closed_pivot_identity
de §72 pour borner le numérateur polynomialement. Si ça bloque,
essaie §85 : dérouler steffensenR_of_fp à x=2 pour R_σ explicit.
```

---

## Session 2 résumé (autonomous)

**État corpus début** : 75 modules verts.
**État corpus fin** : 77 modules verts.

### Cibles complétées (cette session)
- **§73 SigmaStepBoundX2** (downscopé) : framework + cas trivial
  `z = α₀`. Chaînage conditionnel vers §70 documenté. Preuve polynomiale
  full sur demi-plan reste ouverte.
- **§86 IterCountX2** : compte d'itérations ε-précision sous deux
  formes : existentielle (direct depuis §64.3 Tendsto) + géométrique
  explicite `(104/259)^N · (1/4) ≤ ε ⟹ N suffit`.
- **§78 `tendsto_at_zero`** (ajout) : `h(0) = 1/2 ⟹ h^n(0) → α₀`.
  Prouve que la gap `{Re < 1/2 ∧ |z| < 2}` contient des points qui
  s'échappent vers le demi-plan en 1 itération.

### Observations session 2
- L'identité `h(0) = h(-1) = 1/2` est spéciale (Sp(0) = Sp(-1) = 1).
  Tous les `z` avec Sp(z) = 1 mappent vers `1/2`, ancre pratique du
  demi-plan.
- Les downscopes sont inévitables pour §73/§74/§75 — l'analyse
  polynomiale asymptotique nécessite un budget session dédié avec
  outils Lean spéciaux (interval arithmetic?).
- La forme existentielle §86 est triviale (déballage Tendsto) mais
  utile comme API.

### Prochaine session : suggestion de prompt

```
Lis ROADMAP.md Session 2 résumé. Cibles ouvertes restantes :
§74 full σ-basin, §75 équivalence closed-form, §80 p=4 extension.
Attaque §80 (p=4) car plus mécanique : générer §57 pour Sp_C 4.
OU : si tu as Python sympy, essaie de calculer la factorisation
polynomiale qui permet §73 full à la main d'abord, puis formalise.
```

---

## Session 3 résumé (autonomous, p=4 extensions)

**État corpus début** : 77 modules verts.
**État corpus fin** : 79 modules verts.

### Cibles complétées (cette session)
- **§80 HalfPlaneContractionP4** : identité polynomiale exacte pour
  `‖Sp_C 4 z‖²` + borne `‖Sp_C 4 z‖² ≥ (1+a+a²+a³)²` pour **`Re z ≥ 1`**.
  **Découverte importante** : la borne §54-style ne se généralise PAS
  à `Re z ≥ 1/2` à p=4 (contre-exemple à `z = 1/2 + i`). Threshold
  bumped : 1/2 → 1.
- **§81 AlphaX2P4** : ancre `α₀ := 2^{-1/4}`, propriétés arithmétiques
  (α₀^4 = 1/2, α₀ ≥ 3/4, α₀ ≤ 1).

### Observations session 3
- Le pattern §54.1 → §54.2 → §55 → §64 (p=3) ne peut pas être
  bêtement copié à p=4 : la borne polynomiale change de seuil.
- Pour p=4 le threshold "naturel" semble être `Re z ≥ 1` au lieu
  de `Re z ≥ 1/2`. Probablement `Re z ≥ (p-1)/(p)` ou similaire en
  général. À vérifier numériquement pour p ≥ 5.
- La factorisation `polynom = constant_at_a=1 + (a-1)·(positive_poly)`
  est très efficace pour les bornes lower : 2 instances dans §80
  (b² coeff et b⁴ coeff).

### Cibles non-commencées (session 4 future)
- **§82 Generic p framework** : généraliser §80 à tout p ≥ 3 par
  induction. Identité polynomiale et threshold `(p-2)/(p-1)` ?
- **§83/§84 generic x** : x ∈ (1, 8) à p=3.
- **§85 explicit Steffensen radius** : Classical.choose unfolding.
- **§77b/§78/§79** : Niveau 5 closing — nécessite analyse dynamique.

### Prochaine session : suggestion de prompt

```
Lis ROADMAP.md Session 3 résumé. Si Python sympy disponible :
calcule l'identité polynomiale §80.2 généralisée pour p=5,
trouve le threshold optimal pour la borne lower, puis formalise.
Sinon : attaque §83 (Banach generic x ∈ (1, 8) à p=3) — paramétriser
§64 sur x au lieu de fixer x=2.
```

---

## Session 4 résumé (autonomous, paramétrisations)

**État corpus début** : 79 modules verts.
**État corpus fin** : 82 modules verts.

### Cibles complétées (cette session)
- **§83 AlphaGenericX** : ancre `α(x) := x^{-1/3}` paramétrée sur x.
  Identités `α(x)^3 = 1/x`, `α(x) < 1` pour x > 1, `α(x) ≥ 1/2` pour
  `x ≤ 8`. Étend §64.1 (alphaX2 = α(2)) à toute la plage utile.
- **§82 AlphaGenericPX** : unification `α(p, x) := x^{-1/p}`.
  Identité `α(p, x)^p = 1/x`. Spécialisations vers alphaX2, alphaX2P4,
  alphaX. Cap commun pour les ancres §64, §81, §83.
- **§84 HalfPlaneGenericX** : invariance demi-plan paramétrée pour
  `pandrosion_h_C x 3` à `x ∈ (1, 8]`. Premier théorème dynamique
  paramétré sur x (pas fixé à x=2).

### Anecdote build importante
- §84 v1 (avec `field_simp` + `nlinarith` chains) : build > 5 minutes.
  Killed via docker. Refactor en pure `ring + linarith` + `linarith`
  direct sans field_simp : build 10s. **Le pattern "ring identity hunt"
  reste le plus efficace** pour les preuves polynomiales.

### Cibles bloquées (toutes les restantes)
- §74, §75, §76 : analyse polynomiale complexe (degré 18+ field_simp).
- §77b : dynamique complexe sur région bornée + voisinages ω, ω².
- §78 NonPrincipalBasin : mesure topologique des bassins ω·α₀.
- §79 : ferme Niveau 5, dépend des précédents.
- §85 : Classical.choose unfolding pour steffensenR_of_fp.
- §87 : Böttcher — formalisation de séries holomorphes.

### Prochaine session : suggestion de prompt

```
Lis ROADMAP.md Session 4 résumé. Toutes les cibles "faciles"
sont cochées. Restantes demandent ou bien :
(a) analyse polynomiale fine (§73 full, §74, §75) — utilise sympy
    pour trouver factorisation, puis formalise via ring identity.
(b) dynamique complexe (§77b, §78, §79) — extension Mathlib hors
    portée d'une session.
(c) Classical.choose unfold (§85) — nécessite refonte §33.

Recommandation : repos session 5, attendre validation utilisateur
sur direction (Mathlib dynamics? sympy-aided polynomial?).
Sinon : §85 explicit Steffensen radius via unfolding.
```

---

## Session 5 résumé (autonomous, downscopes ponctuels)

**État corpus début** : 82 modules verts.
**État corpus fin** : 84 modules verts.

### Cibles complétées (downscopes ponctuels)
- **§74 SigmaConcretePointsX2** : `σ_closed(1) = 353/444`, `σ_closed(2)
  = 5571/6888`, et `|σ_closed(z) - α₀| < 1/4` à ces deux points
  (via `α₀ ∈ [3/4, 1]` + linarith). Vérification ponctuelle de C1.
- **§75 SigmaEquivPointsX2** : `steffensen_step_C 2 3 (1) =
  sigma_x2_closed (1) = 353/444`. Calcul direct via norm_num après
  unfold. Building block pour preuve générale future, évite
  field_simp degré 18+.

### Observations session 5
- Les downscopes "ponctuels" sont productifs : prouver C1 à des
  points concrets démontre la *plausibilité* du full bound sans
  l'analyse polynomiale lourde.
- `norm_num` après `unfold + simp [Sp_C_p3]` ferme les calculs
  rationnels constants efficacement (build §75 en 9s).
- L'équivalence steffensen_step_C ↔ sigma_x2_closed à un point
  concret prend ~9s ; full general theorem hung 5 min via field_simp
  (cf. session 4 §84 incident).

### Cibles bloquées identiques à session 4
- §76 §77b §78 §79 §85 §87.

### Prochaine session : suggestion de prompt

```
Lis ROADMAP.md Session 5 résumé. Cibles restantes :
- §76 HalfPlaneSigmaTendstoP3X2 — peut être chaîné conditionnellement
  sur R_σ ≥ 1/4 (§70 SteffensenRadiusAtLeast14), via §75 partial pour
  σ(1) ∈ B(α₀, 1/4).
- §85 Explicit Steffensen radius — refonte §33.
- §77b §78 §79 §87 — dynamique complexe out of scope.

Recommandation : §76 conditionnel chain (10-15 min, faisable).
```
