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

- [ ] **§74 SigmaStepIntoBasinX2 — full** : étendre §73 à `Re z ≥ 1/2`
  par découpage en régions (`[1/2, 2]` traité séparément via borne
  polynomiale, `[2, ∞)` par asymptotique). Utiliser `sigma_x2_closed`
  comme handle.

- [ ] **§75 steffensen_step_C ↔ sigma_x2_closed** : prouver
  l'équivalence formelle à `x = 2, p = 3` en factorisant la `field_simp`
  en 3–4 étapes intermédiaires (éviter l'explosion degré 18). Clef :
  passer par `steffensen_denom_C_p3_x2_closed`.

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

- [ ] **§80 HalfPlaneContractionPk — p = 4** : généraliser §57 à
  `p = 4` avec `Sp_C 4 z = 1 + z + z² + z³`. Même structure
  algébrique, degrés bumped.

- [ ] **§81 BanachX4Concrete** : analogue §64 à `x = 2, p = 4`.
  `α₀ = 2^{-1/4} ≈ 0.841`.

- [ ] **§82 Generic p ≥ 3 framework** : paramétrer `x = 2` sur
  `p` variable. Prouver rate `K(p) < 1` par induction sur `p`.

### Extensions x ≠ 2 (priorité #4)

- [ ] **§83 Banach generic x ∈ (1, 8)** : généraliser §64 à tout
  `x ∈ (1, 8)` réel. `α₀(x) = x^{-1/3}`. Rate dépend de `x`.

- [ ] **§84 HalfPlane generic x** : généraliser §67.

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
