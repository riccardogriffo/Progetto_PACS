# LEGGIMI #
PACS project, all rights reserved.

***TODO***

1 - Scrivere una classe problem che applichi tutti gli operators per implementare Laplace

2 - Scrivere pezzi per Dirichlet non omogeneo

3 - Plottare, capire che librerie esistono e che tipo di dati ricevono

4 - Scrivere classe che trasforma il vettore risultato in una valutazione della soluzione su una griglia

    (input : griglia in forma di un doppio vettore su cui valutare la soluzione, soluzione, fespace
    cosigli: utilizzare la valutazione effettuata dalla classe HierarchicalBasis, assegnare a
    ciascun elemento i punti della griglia in ingresso che gli competono
    output: matrice con le valutazioni sul prodotto cartesiano dei vettori in ingresso delle
    ascisse e delle ordinate)

5 - Capire come implementare Neumann e come riciclare il codice scritto per Dirichlet

    (consigli: sostituire al posto della derivata normale di u il valore dato nella condizione 
    al contorno)
6 - Scrivere metodi per calcolare le norme e verificare le stime di convergenza

Per il momento basta così.
