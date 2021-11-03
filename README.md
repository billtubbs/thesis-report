# Package ulthese

The package provides a class based on memoir and compatible with LaTeX
and XeLaTeX to prepare theses and memoirs compliant with the
presentation rules set forth by the Faculty of Graduate Studies of
Université Laval, Québec, Canada. The class also comes with an
extensive set of templates for the various types of theses and memoirs
offered at Laval.

Please note that the documentation for the class and the comments in
the templates are written in French, the language of the target
audience.

## Licence

LaTeX Project Public License, version 1.3c or (at your option) any
later version.

## Version

5.3a (2021-08-09)

## Author

Université Laval <ulthese-dev@bibl.ulaval.ca>

## Source code repository

https://gitlab.com/vigou3/ulthese/

> The rest of this file is in French for the target audience.

# Paquetage ulthese

Le paquetage **ulthese** fournit la classe du même nom permettant de
composer des thèses et des mémoires immédiatement conformes aux règles
générales de présentation matérielle de la Faculté des études
supérieures et postdoctorales (FESP) de l'[Université
Laval](https://www.ulaval.ca). 

La classe permet également de produire les types de documents suivants
selon les mêmes règles de présentation que les thèses et mémoires:
examen de doctorat, essai de maîtrise, projet de recherche, rapport de
stage.

La classe est compatible avec LaTeX et XeLaTeX.

## Contenu du paquetage

- `ulthese.ins`: fichier d'installation de la classe; voir ci-dessous;
- `ulthese.dtx`: fichier source documenté de la classe;
- `ulthese.pdf`: documentation de la classe;
- `ul_p.eps`:    logo de l'Université Laval en format EPS;
- `ul_p.pdf`:    logo de l'Université Laval en format PDF;
- `README.md`:   le présent fichier.

## Installation

Pour installer la classe ulthese, compiler avec LaTeX le fichier
`ulthese.ins`, en exécutant

    latex ulthese.ins

depuis une ligne de commande. On peut aussi ouvrir le fichier
`ulthese.ins` dans son éditeur LaTeX et lancer la compilation avec
LaTeX, pdfTeX, pdfLaTeX ou XeLaTeX depuis celui-ci.

La compilation du fichier ulthese.ins générera plusieurs fichiers:

1. la classe elle-même:
   - `ulthese.cls`

2. des gabarits pour le document maître de différents types de
   thèses et de mémoires:
   - `gabarit-doctorat.tex`
   - `gabarit-doctorat-articles.tex`
   - `gabarit-doctorat-cotutelle.tex`
   - `gabarit-maitrise.tex`
   - `gabarit-maitrise-articles.tex`
   - `gabarit-maitrise-bidiplomation.tex`

3. des gabarits pour quelques parties d'un document:
   - `resume.tex`
   - `abstract.tex`
   - `remerciements.tex`
   - `avantpropos.tex`
   - `introduction.tex`
   - `chapitre1.tex`
   - `chapitre1-articles.tex`
   - `chapitre2.tex`
   - `chapitre2-articles.tex`
   - `conclusion.tex`
   - `annexe.tex`

Les gabarits inutiles peuvent être supprimés.

## Documentation

Le fichier `ulthese.pdf` contient la documentation complète de la
classe et des gabarits. Le document peut être recréé à partir du code
source avec les commandes

    pdflatex ulthese.dtx
    makeindex -s gglo.ist -o ulthese.gls ulthese.glo
    pdflatex ulthese.dtx

## Historique des versions

L'historique des versions de la classe se trouve en annexe de la documentation.

## Commentaires et suggestions

Le [dépôt du projet](https://gitlab.com/vigou3/ulthese/) dans GitLab
demeure le meilleur endroit pour rapporter des bogues ou pour proposer
des améliorations à la classe.

## Aide additionnelle

Faire appel aux ressources suivantes, dans l'ordre:

1. consulter le [WikiThèse de l'Université Laval](http://www.theses.ulaval.ca/wiki/);

2. consulter les archives de la liste de distribution
   [ulthese-aide](http://listes.ulaval.ca/listserv/archives/ulthese-aide.html);

3. écrire à la liste de distribution [ulthese-aide](mailto:ulthese-aide@listes.ulaval.ca).
