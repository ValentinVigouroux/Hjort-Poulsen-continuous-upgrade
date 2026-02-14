clear all 
set more off 

cd "C:\Users\valen\Documents\MAG3\EVALUATION IMPACT Lopers\TD\"

use "grid_panel_final.dta" 

keep if dist_backb <= 20 // on exclue toute zone éloignée de plus de 20 km du back bone dans l'étude.  "  les contrôles sont les cellules situées entre 10 et 20 km" la phrase restreintes à l'échantillon in_sample ==1 laisse entendre que toutes les observations à < 20 km du backbone ont in_sample ==1 Mais ce n'est pas le cas. 

preserve 
    collapse (mean) night, by(year)
twoway line night year, sort ///
        xlabel(2000(1)2019, angle(45) grid) ///
        title("Évolution de la moyenne de l'intensité lumineuse (2000-2019)") ///
        xtitle("Année") ytitle("Moyenne de Night") 
restore 
 


cap drop country_year // saute cette étape si déjà généré 
egen country_year = group(country year)
// i.country_year = effets fixes pays années

// j'utilise la méthode utilisée dans le papier de 2019 et vue en cours avec areg 

// absorb() permet d'inclure des effets fixes lorsqu'il y a un grand nombre de "catégories"  effet fixe cellule : 88K cellules = 88 k variables = beaucoup 

//contrairement au auteurs, nous n'allons pas différencier les cellules traitées (<10km) selon leur distance du backbone 

areg ihs_night_dmsp did i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)
lincom did  // permet de n'afficher que la première ligne du tableau 


areg ln_gdp did i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)
lincom did 

areg ihs_gdp did i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)
lincom did

label variable ihs_night_dmsp "Night DMSP" //je renomme l'étiquette de l'a variable car c'était trop long
eststo clear 

eststo r1: areg ihs_night_dmsp did access pop_final trend_proxy i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)

eststo r2: areg ln_gdp did access pop_final trend_proxy i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)

eststo r3: areg ihs_gdp did access pop_final trend_proxy i.country_year if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)

//LATEX 
esttab r1 r2 r3 using "resultats_did.tex", replace ///
    label b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(did access pop_final trend_proxy) /// retirer les lignes country_year
    stats(N r2_a, fmt(%12.0f %12.4f) labels("Observations" "R2-ajusté")) ///
    title("Effet de l'arrivée des câbles sous-marins") ///
    booktabs nonotes

	/*
// Exportation vers word (facultatif automatisation)
esttab r1 r2 r3 using "resultats_did.rtf", replace ///
    label /// 
    b(4) se(4) /// 
    star(* 0.10 ** 0.05 *** 0.01) /// Seuils de significativité standard 
    keep(did access pop_final trend_proxy) /// Masque les effets fixes pays-année 
    order(did access pop_final trend_proxy) /// Ordonne les variables proprement
    stats(N r2_a, labels("Observations" "R2-ajusté")) /// Statistiques de bas de tableau 
    title("Tableau 1: Effet de l'arrivée des câbles sous-marins") /// Titre du tableau 
    modelwidth(10) /// Fixe une largeur de colonne pour éviter les décalages dans Word
    varwidth(20) /// Largeur pour le nom des variables
    nogaps /// Supprime les espaces inutiles entre les lignes
    onecell /// Met le coefficient et l'erreur standard dans la même cellule (optionnel, plus compact)
    addnotes("Note : Effets fixes cellule et pays-année inclus. Erreurs-types clusterisées au niveau cellule.")
*/ 	


// on restreint l'échantillon à l'intérieur de la commande pour la suite 


estpost tabstat base_night ln_elev ln_base_pop if year == 2000, by(is_close) statistics(mean sd) columns(statistics)

// Exportation
esttab . using "Moyenne_1.tex", replace ///
    cells("mean(fmt(3)) sd(par fmt(2))") /// 
    label ///                                 
    booktabs ///                              
    unstack ///                               
    nonumber ///                              
    mtitles          // 
	
* Test Balance 
// Pour tester une différence entre les échantillons contrôle et traitement  on peut recourir à la regression suivante : 
eststo clear
eststo: regress base_night is_close if year == 2000
eststo: regress ln_elev is_close if year == 2000
eststo: regress ln_base_pop is_close if year == 2000

esttab using "balance1.tex", replace ///
    label /// 
    b(3) se(3) ///  
    star(* 0.10 ** 0.05 *** 0.01) /// 
    nodepvars ///            
    mtitles("Luminosité" "Altitude" "Population") /// 
    stats(N, fmt(%12.0f %12.4f) labels("Observations")) ///    
    booktabs                                  
	
      // PSM 
	  
// l'échantillon est restreint car nous souhaitons matcher les individus avant le traitement pour ne pas "contaminer" notre matching 

probit is_close base_night ln_elev ln_base_pop ln_base_gdp is_agri is_forest  ln_access if year == 2000 

* 2. Création du PSCORE
cap drop PSCORE
predict PSCORE, pr

/*

// Test du Common Support (facultatif)
twoway (kdensity PSCORE if is_close==1, lcolor(blue) lwidth(medthick)) ///
       (kdensity PSCORE if is_close==0, lcolor(red) lpattern(dash)), ///
       legend(order(1 "Traités" 2 "Non traités")) ///
       title("Densité du score de propension") ///
       xtitle("Score de propension") ytitle("Densité")
*/ 


//  MATCHING  j'utilise à présent psmatch2  

// Cette étude est caractérisée par un très grand nombre d'observations N = 88000 . 
//Je choisis d'utiliser un Caliper/radius matching avec une bande faible et un voisin le plus proche, cette restriction à un seul voisin très proche en p-score minimise le biais.

// Nous disposons d'un grand nombre d'observations donc nous pouvons nous permettre de "jeter" certaines observations qui serait trop distantes en p-score
// De même, nous effectuons un matching avec remise càd que Stata peut réutiliser plusieurs fois la meme cellule controle pour réaliser le matching. 
// Cela permet de minimiser d'avantage le biais au détriment de la puissance. 

psmatch2 is_close if year == 2000, pscore(PSCORE) outcome(ln_gdp) neighbor(1) caliper(0.001) common 

cap drop support_cell // necessaire pour relancer le code. écrase la variable si déjà générée 
cap drop weight_cell
bysort cell_id: egen support_cell = max(_support) // pour chaque cellule si _support prend la valeur 1 durant une année elle va générer 1 pour support _cell => étendre les données pour les cellules matchés à toutes les années 
bysort cell_id: egen weight_cell = max(_weight) // nous étendons également les données de poids qui prennent en comptes les observations controles matchés plusieurs fois. (matching avec remise)


keep if support_cell == 1 // on ne garde que les individus matchés (années 2000 et +)
eststo clear

// did cell= is_close =1 ET n_cables > 1 
eststo psm_did: quietly areg ln_gdp did access pop_final trend_proxy i.country_year [pweight=weight_cell] if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)


esttab psm_did using "Resultats_PSM_DiD.tex", replace ///
    label /// 
    b(4) se(4) /// 
    star(* 0.10 ** 0.05 *** 0.01) ///
    keep(did access pop_final trend_proxy) /// 
    order(did access pop_final trend_proxy) /// 
    stats(N r2_a, fmt(%12.0f %12.4f) labels("Observations" "R2-ajusté")) /// 
    title("Impact de l'Internet rapide sur le PIB (Estimation PSM-DiD)") ///
	mtitles("ln(gdp)") ///
    booktabs /// 
    nonotes 
   
// version rtf 
/*
	esttab psm_did using "Resultats_PSM_DiD.rtf", replace ///
    label /// 
    b(4) se(4) /// 
    star(* 0.10 ** 0.05 *** 0.01) ///
    keep(did access pop_final trend_proxy) /// Masque les variables pays années
    order(did access pop_final trend_proxy) ///
    stats(N r2_a, fmt(%12.0f %12.4f) labels("Observations" "R2-ajusté")) ///
    title("Tableau 2 : Impact de l'Internet rapide sur le PIB (Estimation PSM-Did)") ///
    mtitles("ln(gdp)") /// 
    varwidth(25) /// 
    modelwidth(12) /// 
    nogaps /// 


	*/ 
	// Spécification améliorée : Intensité de traitement 	

gen did_intensity= post*n_cables
eststo clear

eststo intensity: areg ln_gdp did did_intensity access pop_final trend_proxy i.country_year [pweight=weight_cell] if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)


esttab intensity using "Resultats_Intensite.tex", replace ///
    label /// 
    b(4) se(4) /// 
    star(* 0.10 ** 0.05 *** 0.01) ///
    keep(did did_intensity access pop_final trend_proxy) /// 
    order(did did_intensity access pop_final trend_proxy) /// 
    stats(N r2_a, fmt(%12.0f %12.4f) labels("Observations" "R2-ajusté")) /// 
    title("Table 3 Impact de l'intensité de l'internet rapide sur le PIB (PSM-Did)") ///
    mtitles("ln(gdp)") ///
    booktabs nonotes

* 3. Exportation vers Word (RTF)
esttab intensity using "Resultats_Intensite.rtf", replace ///
    label /// 
    b(4) se(4) /// 
    star(* 0.10 ** 0.05 *** 0.01) ///
    keep(did did_intensity access pop_final trend_proxy) /// 
    order(did did_intensity access pop_final trend_proxy) /// 
    stats(N r2_a, labels("Observations" "R2-ajusté")) /// 
    title("Table 3 Impact de l'intensité de l'internet rapide sur le PIB (PSM-Did)") ///
    mtitles("ln(gdp)") /// 
    varwidth(30) /// 
    modelwidth(12) /// 
    nogaps /// 
    addnotes("Note : (PSM). did_intensity capte l'effet marginal d'un cable supplémentaire")	
	
// Version 2 : intensité 

gen d_1_cable = (n_cables == 1 & post == 1 & is_close == 1)
gen d_2_plus   = (n_cables >= 2 & post == 1 & is_close == 1)

/*
* Régression intensité2 
eststo intensity_robust: areg ln_gdp d_1_cable d_2_plus access pop_final trend_proxy i.country_year [pweight=weight_cell] if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)
*/ 
	
// Test placebo 

// Tout d'abord j'ai vérifié qu'il y avait une date unique de premier cable sous-marin pour chaque pays en faisant :  tab country cable_year 

cap program drop test_placebo // on nettoie les reste du précédent programme pour le relancer 
program test_placebo, rclass // program permet de définir une suite d'instructions que je vais executer x fois avec simulate 
    preserve
        * 1. On prépare les fausses dates au niveau pays
        collapse (first) cable_year, by(country) // on garde seulement les dates de 1er raccordement sous marin par pays 
		
		// séquence de commande pour mélanger les années => abandonnée car les années sont trop proche forte correlation qui limite l'efficacité de notre test placebo 
		
/*		gen r = runiform() // genere une valeur aléatoire dans une colonne
		putmata x = cable_year, replace // copie la commande cable year dans une matrice appelée x 
		sort r  // trie le dataset selon r (donc aléatoirement)
		getmata rand_cy = x, replace // prend les valeur stockées (avant le mélange dans la matrice et les injecte dans le dataset via la variable rand_cy
		*/ 
		
		gen rand_cy = runiformint(2000, 2013)
        keep country rand_cy
        tempfile placebo_d
        save `placebo_d'
    restore

    merge m:1 country using `placebo_d', nogenerate

    gen P_post = (year >= rand_cy)
    gen P_did = P_post * is_close
    
    // regression 
	keep if support_cell == 1 // modèle PSM DiD
    quietly areg ln_gdp P_did access pop_final trend_proxy i.country_year [pweight=weight_cell] if year >= 2000 & year <= 2013, a(cell_id) cluster(cell_id)
    
    // On "extrait" le coefficient pour que simulate puisse le récupérer
    return scalar c_placebo = _b[P_did]
    
    //nettoyage des variables pour la prochaine itération
    drop rand_cy P_post P_did
end

// lancer la boucle qui reprend le programme ci dessus 
simulate c_placebo = r(c_placebo), reps(50) seed(12345): test_placebo
histogram c_placebo, title("Test Placebo (50 itérations)") ///
    xline(0.002, lcolor(red)) // 
	
   
	   
	
