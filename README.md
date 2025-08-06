Here you can find all the Python scripts that I have used for my master's thesis 2025. The link below will take you to the university's database where you can find the full paper:

https://uu.diva-portal.org/smash/record.jsf?dswid=6492&pid=diva2%3A1975073&c=2&searchType=SIMPLE&language=sv&query=Simon+Liljeblad&af=%5B%5D&aq=%5B%5B%5D%5D&aq2=%5B%5B%5D%5D&aqe=%5B%5D&noOfRows=50&sortOrder=author_sort_asc&sortOrder2=title_sort_asc&onlyFullText=false&sf=all

Important note: 
ChatGPT HAS been used as a supporting tool for creating some of the code. This, however, is within the boundaries of what is allowed by Uppsala University (as of the making of the report).

Below follows a quick explanation of what every script does

PCA.ipynb: 
This is the main script used to perform dimensional reduction using PCA

TSNE.ipynb: 
Performs dimensional reduction using TSNE

UMAP.ipynb: 
Performs dimensional reduction using UMAP

aligning.ipynb: 
This script aligns the unaligned structures to the first-frame wild-type

bars.ipynb: 
This script takes the maximum ARI score and creates nice bar-charts including the used dimensional reduction method and photon energy. These bar charts are then fused together in an external editing software. You can see my report if you want to know how they look like.
(Note here that the scores are input by the user directly, which I did by analyzing the score images in /home/simon/3d_plots/)

charges.ipynb: 
This script calculates the mean charges of every atom-type in the ubi-mutants for every photon energy used. These are presented as a bar-chart.

cross_section.ipynb: 
Used to create plot containing the cross-section of all atom-types in the ubi-mutants, dependant on the photon energy. This also includes the weighted cross section for the average ubi-mutant. You can read more about this in the report, which includes what sources were used.

detector.ipynb: 
The main script used for observing explosion patterns from every single simulation, including an averaged pattern between all simulation.

energies.ipynb: 
Plots the kinetic energies over the entire simulation time for every photon energy used.

gaussian_and_Morse.ipynb: 
This was simply the script used to make the Gaussian-function and the Morse-function used in the report.

mean_charge.ipynb: 
Plots the mean charge over the entire simulation time for every photon energy simulated. Also plots the pulse to get a better understanding of how the increase in charges relate to the pulse.

run_simulations.py: 
The MAIN script used to run simulations, given certain parameters inputed in the script. Why is this file a .py-file and not a .ipynb-file? It's simply so that I could run it with nohup (the simulations can take many hours).

Thank you for reading this. I hope you have a wonderful day :)
If you have any question, you can contact me via mail: liljeblad.simon@gmail.com
