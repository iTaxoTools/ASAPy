/*
	Copyright (C) 2015-2016 G Achaz/ S Brouillet

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

for more information, please contact guillaume achaz <guillaume.achaz@mnhn.fr>/<sophie.brouillet@mnhn.fr>

*/
/******
        file     : asap_common
        function : all fns neened by asap (web version and command line version) exept the core algorithm


        created  : November 2015


        author   : madamesophie


*****/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <ctype.h>
#include <limits.h>
#include "asap.h"
#ifdef ASAP_CL
#define WORKDIR ""
#endif
//#include "drawMat.h"
#include "gdtosvg.h"

#define NBTEST 10000
#ifdef MACOSX
#include <float.h>
#elif defined(SGI)
#include <limits.h>
#elif defined(LINUX)
#include <values.h>
#endif

#ifdef _PROTOTYPES
#undef _PROTOTYPES
#include <float.h>
#define _PROTOTYPES
#else
#include <float.h>
#endif

unsigned long idum_ran1 = 0;
char *strcasestr(const char *haystack, const char *needle);


char debug;


/****************************************************/
/*--------------------------------------------------*/
/* Misc functions*/
/*--------------------------------------------------*/
void swap(int *t1, int *t2)
{
	int v;
	v = *t1;
	*t1 = *t2;
	*t2 = v;
}
void exit_properly(char *ledir)
{
	char commande [1024];
	if (strlen (ledir) >1)
	{
	sprintf (commande, "mv %s%s/results_.html %s%s/results.html", WORKDIR, ledir, WORKDIR, ledir);

	system(commande);
	}
}
/*--------------------------------------------------*/
/*COMPARAISON FNS*/
/*--------------------------------------------------*/
/*ascending sort on distance*/
int compareCase(void const *a, void const *b)
{
	DistPair *pa = (DistPair  *) a;
	DistPair *pb = (DistPair  *) b;

	double x = pa->d;
	double y = pb->d;

	if (x < y)
	return (-1); 
	else if (x > y)
		return (1);
	else

	{
		if (pa->a == pb->a)
			return (pa->b - pb->b) ;
		else
			return (pa->a - pb->a) ; 
	}

}
/*--------------------------------------------------*/
/*descending sort on slope*/
int compareParameter(void const *a, void const *b)
{
	Results *pa = (Results  *) a;
	Results *pb = (Results  *) b;

	double x = pa->other_parameter;
	double y = pb->other_parameter;
//printf("%lf %lf  \n",x,y);
	if (x > y) return (-1);
	else return (1);
}

/*--------------------------------------------------*/
/*ascending sort on probabiliy*/
int compareProba(void const *a, void const *b)
{
	Results *pa = (Results  *) a;
	Results *pb = (Results  *) b;

	double x = pa->proba;
	double y = pb->proba;

	if (x < y) return (-1);
	else return (1);
}

/*--------------------------------------------------*/
/*ascending sort on asap rank*/
int compareRang(void const *a, void const *b)
{
	Results *pa = (Results  *) a;
	Results *pb = (Results  *) b;

	double x = pa->score;
	double y = pb->score;


	double x2=pa->proba;
	double y2=pb->proba;


	double x1=pa->other_parameter;
	double y1=pb->other_parameter;
//printf("%lf %lf  \n",x,y);
	if (x < y) 
		return (-1);
	else 
	{
		if (x > y)
		 return (1);
		
		else //egalite des scores on teste les proba
			{
			if (x2 < y2) return (-1);
				else 
					if (x2 >y2)return (1);
				else //egalité des probas
					if (x1 > y1) return (-1);
			}
	}
		return (1);
}


/*--------------------------------------------------*/
/*descending sort on asap score*/
int compareScore(void const *a, void const *b)/*Plus grand au plus petit*/
{
	Results *pa = (Results  *) a;
	Results *pb = (Results  *) b;	

	double x = pa->score;
	double y = pb->score;
//printf("%lf %lf  \n",x,y);
	if (x > y) return (-1);
	else return (1);
}

/*--------------------------------------------------*/
/*descending sort on distance*/
int compareSpecies(void const *a, void const *b)/*Plus grand au plus petit*/
{
	Results *pa = (Results  *) a;
	Results *pb = (Results  *) b;	

	double x = pa->nbgroups;
	double y = pb->nbgroups;
//printf("%lf %lf  \n",x,y);
	if (x > y) return (-1);
	else return (1);
}

/*ascending sort on graphical position*/
int compareCoord(void const *a, void const *b)
{
	int *pa = (int  *) a;
	int *pb = (int  *) b;

	double x = *pa;
	double y = *pb;
//printf("%lf %lf  \n",x,y);
	if (x < y)

		return (-1); //on trie sur distance
	else 
		return (1);

	

}



/****************************************************/
/*--------------------------------------------------*/
/* Some math function used later*/
/*--------------------------------------------------*/
// [0,1[
double unirandom(){

	int U=RAND_MAX;
	while(U==RAND_MAX)
		U=rand();
		
	return U/(double)RAND_MAX; 
}

/*--------------------------------------------------*/
// [min,max], tinily in favor of min
int uniInt(int min,int max){
	return (rand()%(max-min+1))+min;	
}

/*--------------------------------------------------*/
double exponentialdev()
{
	double dum;

	do
		dum = unirandom();
	while ( dum == 0.0 );
	return -log(dum);
}



/*--------------------------------------------------*/

double poissondev(double mean)
{
	int em;
	double g, t;

	g = exp(-mean);
	em = -1;
	t = 1.0;
	do {
		em += 1;
		t *= unirandom();
	} while (t > g);

	return ((double)em);
}


/*--------------------------------------------------*/
//PRINT FNS ONLY FOR DEBUGGING // VERIF
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/******************************/
void PrintSize( LeftRight *size, int n, int min, int max ){
	
	int j;
	
	for(j=0;j<n;j++)
		printf("%2d,%2d   ",size[j].L,size[j].R);
	printf(" ; min: %d max: %d\n", min, max);
}
void PrintLength( double *lengthTreeLeft,double *lengthTreeRight, int n )
{
	int i;

	for (i=0;i<n;i++)
		printf("%f, ",lengthTreeLeft[i]);
	printf("\n");

	for (i=0;i<n;i++)
		printf("%f, ",lengthTreeRight[i]);
	printf("\n");
}
void print_comp(Composante cc, int nb, Tabcompo *strucompo)
{
	int i, j;
	printf("nb de composante :%d \n", cc.nc);


	for (i = 0; i < nb; i++)
	{

		printf("comp %d (id=%d) (nbelts=%d) :", i, cc.node_compid[i], cc.n_in_comp[i]);

		for (j = 0; j < cc.n_in_comp[i]; j++)
			printf("%d \t", cc.comp[i][j]);
		printf("\nstrucompo[%d].nb=%d:", i, strucompo[i].nb);
		for (j = 0; j < strucompo[i].nb; j++)
			printf("%d\t", strucompo[i].effcompo[j]);
		printf("\n");
	}
	printf("\n***********************\n");
}


/*--------------------------------------------------*/
void print_compNames(Composante cc, DistMat mat)
{
	int i, j;
	int nb = mat.n;
	printf("nb de composante :%d \n", cc.nc);


	for (i = 0; i < nb; i++)
	{
		if (cc.n_in_comp[i] > 1)
		{
			printf("comp %d (id=%d) (nbelts=%d) :", i, cc.node_compid[i], cc.n_in_comp[i]);

			for (j = 0; j < cc.n_in_comp[i]; j++)
				printf("%s \t", mat.names[cc.comp[i][j]]);
			printf("\n");
		}
	}
	printf("\n***********************\n");
}

/*--------------------------------------------------*/
void fprint_htmlcomp(Composante cc, int nb, FILE *f, int *nno)
{
	int i, j;
	fprintf(f, "nb de composante :%d <BR>\n", cc.nc);


	for (i = 0; i < nb; i++)
	{
//		double nn= (cc.n_in_comp[i]*(cc.n_in_comp[i]-1))/2.0;
		if (cc.n_in_comp[i] >= 1)
		{
			fprintf(f, "comp %d (id=%d) (nbelts=%d) :", i, cc.node_compid[i], cc.n_in_comp[i]);

			for (j = 0; j < cc.n_in_comp[i]; j++)
				fprintf(f, "%d \t", cc.comp[i][j]);
			fprintf(f, "<BR>\n");
			if (nno != NULL)
				for (j = 0; j < cc.n_in_comp[i]; j++)
					fprintf(f, "%d \t", nno[cc.comp[i][j] ]);
			fprintf(f, "<BR>\n");
		}
	}
	fprintf(f, "<BR>\n***********************\n<BR>");
}

/*--------------------------------------------------*/
void print_comp_file(Composante cc, int nb, FILE *f)
{
	int i, j;
	fprintf(f, "nb de composante :%d \n", cc.nc);


	for (i = 0; i < nb; i++)
	{
//		double nn= (cc.n_in_comp[i]*(cc.n_in_comp[i]-1))/2.0;

		fprintf(f, "comp %d (id=%d) (nbelts=%d) :", i, cc.node_compid[i], cc.n_in_comp[i]);

		for (j = 0; j < cc.n_in_comp[i]; j++)
			fprintf(f, "%d \t", cc.comp[i][j]);

		fprintf(f, "\n");
	}
	fprintf(f, "\n***********************\n");
}


/*--------------------------------------------------*/
void print_tab(int *t, int nb)
{
	int i;
	printf("%d elts:", nb);
	for (i = 0; i < nb - 1; i++)
		printf("%d,", t[i]);
	printf("%d\n", t[i]);

}



/*--------------------------------------------------*/
int compte_comp(Composante cc, int nb)
{
	int i, nn = 0;
	printf("nb de composantes :%d \n", cc.nc);


	for (i = 0; i < nb; i++)
	{
//		double nn= (cc.n_in_comp[i]*(cc.n_in_comp[i]-1))/2.0;
		if (cc.n_in_comp[i] != 0)
			nn++;
	}
	if (cc.nc != nn)
		printf("\n**********aaararraragazrazraztgargzar ***********************\n");
	return (nn);
}




/****************************************************/
/*--------------------------------------------------*/
/*INIT , REINIT OR CLEAN FONCTIONS*/
/*--------------------------------------------------*/
void reinit_nod(int nleaves, Node *simnodes)
{
	int i;
	int nbnodes = (2 * nleaves) - 1;
	for (i = 0; i < nleaves; i++)
	{

		simnodes[i].nbdesc = 0;
		simnodes[i].anc = -1;
		simnodes[i].time = 0;
		simnodes[i].nbmut = 0;
		simnodes[i].nb_under = 1;
		
	}

	for (i = nleaves; i < nbnodes; i++)
	{
		simnodes[i].nbdesc = 0;
		simnodes[i].anc = -1;
		simnodes[i].time = 0;
		simnodes[i].nbmut = 0;
		simnodes[i].nb_under = 0;
	}
}


/*--------------------------------------------------*/
/*init struct tabcompo*/
/*
	Now it is a [nbseq * nbseq] matrix
	filled with -1, with the exception of the first cell that is 1 because init state is one compo

*/
void inittabcompo(Tabcompo *strucompo, int nbseq, FILE *ff, char *ledir)
{
	int i;
	if (ff == NULL)
		ff = stderr;
	for (i = 0; i < nbseq; i++)
	{
		strucompo[i].nb = 1;
		strucompo[i].effcompo = (int *)malloc(sizeof(int) * nbseq);
		if (strucompo[i].effcompo == NULL) {
			fprintf(ff, "inittabcompo: MEMORY ERROR cant alloc strucompo[%d].effcompo\n", i);
			if (ff != NULL)
				{fclose(ff); exit_properly(ledir);}
			else exit(1);
		}
		memset(strucompo[i].effcompo, -1, nbseq);
		strucompo[i].effcompo[0] = 1;
		strucompo[i].nodecompo = (int *)malloc(sizeof(int) * nbseq);
		if (strucompo[i].nodecompo == NULL) {
			fprintf(ff, "inittabcompo: MEMORY ERROR cant alloc strucompo[%d].effcompo\n", i);
			if (ff != NULL)
				{fclose(ff); exit_properly(ledir);}
			else exit(1);
		}
		memset(strucompo[i].nodecompo, -1, nbseq);
		strucompo[i].nodecompo[0] = i;

	}
}

/*--------------------------------------------------*/
/*
	
*/
void initNodes(FILE *f, Node *zenodes, struct DistanceMatrix mat, char *ledir)
{
	int i;
	int nbnodes = (2*mat.n) - 1 ;

	for (i = 0; i < nbnodes; i++)
	{

		zenodes[i].anc = -1;

		zenodes[i].x = -1;

		zenodes[i].y = -1;

		zenodes[i].dist = 0;

		zenodes[i].name = NULL;
		zenodes[i].pval = 0;

		zenodes[i].nbdesc = 0;
		zenodes[i].sum_inter = 0;
		zenodes[i].sum_all = 0;
		zenodes[i].desc = NULL;
		zenodes[i].nb_inter = 0;
		zenodes[i].nb_all = 0;
		//zenodes[i].tmrca_obs = 0;
		zenodes[i].round = -1;
		zenodes[i].nb_under = 0;
		zenodes[i].color=(i+1)%16;//assign a color 
		zenodes[i].first_to_draw=i;
		zenodes[i].to_be_checked=1;

		zenodes[i].S_all_theo=0.0;
		zenodes[i].S_intra_theo=0.0;

	}


	for (i = 0; i < mat.n; i++)
	{
		zenodes[i].name = malloc((strlen(mat.names[i]) + 1) * sizeof(char));
		if (zenodes[i].name == NULL)
		{	fprintf(f, "ERRORRRR<BR>\n"), fclose(f), exit_properly(ledir);
			if (f != NULL)
			{fclose(f); exit_properly(ledir);}
			else exit(1);
		}

		strcpy(zenodes[i].name, mat.names[i]);
		zenodes[i].color=i%16;
		zenodes[i].nb_under = 1;

	}

}
/*--------------------------------------------------*/
/*clear  all struccompo which have been modified */
void clearalltab(Tabcompo *strucompo, Composante *comp, int nbseq)
{
	int i, j;


	for (j = 0; j < comp->naltered; j++)
	{
		i = comp->altered_comp[j];
		if (i == -1) fprintf(stderr, "this is the end of the world: negative value saved\n"), exit(1);
		strucompo[i].nb = comp->n_in_comp[i];
		memset(strucompo[i].effcompo, -1, nbseq);
		strucompo[i].effcompo[0] = comp->n_in_comp[i];
		strucompo[i].nb = 1;
		comp->altered[i] = 0;
	}
	memset(comp->altered_comp, -1, comp->naltered);
	memset(comp->altered,0,nbseq);
	comp->naltered = 0;
}


/*--------------------------------------------------*/
void cleartabcompo(Tabcompo *strucompo, Composante comp, int nbseq)
{
	int i;

	for (i = 0; i < nbseq; i++)
		if (strucompo[i].nb != 1)
		{
			strucompo[i].nb = comp.n_in_comp[i];
			memset(strucompo[i].effcompo, -1, nbseq);
			strucompo[i].effcompo[0] = comp.n_in_comp[i];
			strucompo[i].nb = 1;
			comp.altered[i] = 0;

		}

	memset(comp.altered_comp, -1, comp.naltered);
	comp.naltered = 0;


}

/*--------------------------------------------------*/
void freecomp(Composante *comp, int nbseq)
{
	int i;


	free(comp->node_compid);

	free(comp->n_in_comp);
	
	free(comp->Sall_in_comp);

	free(comp->altered_comp);

	free(comp->altered);

	comp->naltered = 0;


	for (i = 0; i < nbseq; i++)
		if(comp->altered_comp_prev[i] != NULL)
			free(comp->altered_comp_prev[i]);

	free(comp->altered_comp_prev);
		
	for (i = 0; i < nbseq; i++)
		free(comp->comp[i]);

	free(comp->comp);
}


/*--------------------------------------------------*/
/* as the name says: init everything in the struc Composante*/
void initcomp(Composante *comp, int nbseq, FILE *ff, char *ledir)
{
	int i, j;
	if (ff == NULL)
		ff = stderr;

	comp->nc = nbseq;         /* number of composantes. At the start each sequence is in a different composante */

	comp->node_compid = (int *)malloc(sizeof(int) * nbseq);
	if (!comp->node_compid) fprintf(ff, "initcomp: MEMORY ERROR error can allocate node_compid\n"), fclose(ff), exit_properly(ledir);

	comp->n_in_comp = (int *)malloc(sizeof(int) * nbseq);
	if (!comp->n_in_comp) fprintf(ff, "initcomp:MEMORY ERROR error can allocate n_in_comp\n"), fclose(ff), exit_properly(ledir);

	comp->Sall_in_comp = (double *)malloc(sizeof(double) * nbseq);
	if (!comp->Sall_in_comp) fprintf(ff, "initcomp:MEMORY ERROR error can allocate Sall_in_comp\n"), fclose(ff), exit_properly(ledir);

	comp->comp = (int **)malloc(sizeof(int*)*nbseq);
	if (!comp->comp) fprintf(ff, "initcomp:MEMORY ERROR error can allocate n_incomp\n"), fclose(ff), exit_properly(ledir);

	comp->naltered = 0;
	comp->altered_comp = (int *)malloc(sizeof(int)*nbseq);
	if (!comp->altered_comp) fprintf(ff, "initcomp:MEMORY ERROR error can allocate n_inaltered_comp\n"), fclose(ff), exit_properly(ledir);

	comp->altered_comp_prev=(int **)malloc(sizeof(int *)*nbseq);
	if (!comp->altered_comp_prev) fprintf(ff, "initcomp:MEMORY ERROR error can allocate altered\n"), fclose(ff), exit_properly(ledir);


	comp->altered = (char *)malloc(sizeof(char) * nbseq);
	if (!comp->altered) fprintf(ff, "initcomp:MEMORY ERROR error can allocate altered\n"), fclose(ff), exit_properly(ledir);



	for (i = 0; i < nbseq; i++)
	{
		comp->comp[i] = (int *)malloc(sizeof(int) * nbseq);
		if (!comp->comp[i]) fprintf(ff, "initcomp:MEMORY ERROR error can allocate n_incomp[%d]\n", i), fclose(ff), exit_properly(ledir);
		for (j = 1; j < nbseq; j++)
			comp->comp[i][j] = -1;
		comp->altered_comp[i] = -1;
		comp->comp[i][0] = i;          /* at the begining, each composante has only its own sequence. They share the same id */
		comp->n_in_comp[i] = 1;
		comp->Sall_in_comp[i] = 0.0;
		comp->node_compid[i] = i;
		comp->altered[i] = 0;
		
		comp->altered_comp_prev[i]=NULL;


	}
}
/*--------------------------------------------------*/
void resetcomp(Composante *comp, int nbseq)
{
	int i, j;

	comp->nc = nbseq;         /* number of composantes */


	comp->naltered = 0;

	for (i = 0; i < nbseq; i++)
	{

		for (j = 1; j < nbseq; j++)
			comp->comp[i][j] = -1;
		comp->altered_comp[i] = -1;
		comp->comp[i][0] = i;
		comp->n_in_comp[i] = 1;
		comp->Sall_in_comp[i] = 0.0;
		comp->node_compid[i] = i;
		comp->altered[i] = 0;
		
		if(comp->altered_comp_prev[i] != NULL)
			free(comp->altered_comp_prev[i]);

		comp->altered_comp_prev[i]=NULL;
	}
}



/*--------------------------------------------------*/
/*
Sum of all nodes excepted the new ones
*/
void sum_intra_notused(int nb, int r, Node *zenodes, double *s, long *nbi)
{
	int i;

	for (i = 0; i < nb; i++)
		if (zenodes[i].round != r && zenodes[i].round != -1)
		{
			*s = *s + zenodes[i].sum_intra;
			*nbi = *nbi + zenodes[i].nb_intra;
		}

}
/*--------------------------------------------------*/
/* Get which desc has the bigger group size under it*/
int get_descmax(Node *n, int x)

{
int i,k,j,m;

i=n[x].desc[0];// le premier desc
m=n[i].nb_under; //initialise m avec le nbr dessous premier desc
j=0;
for (i=1;i<n[x].nbdesc;i++)
{
	k=	n[x].desc[i];
	if (m<n[k].nb_under)
		{m=n[k].nb_under;
		j=i;}
}	
return j;	

}

/*--------------------------------------------------*/
/* associate a color to a node, usefull for nice drawings*/
void color_clado(Node *zenodes, int current_node ,int *whichcolor)

{
int i;
int n;
//int color;

	zenodes[current_node].color= *whichcolor;

	if(zenodes[current_node].nbdesc==0)
		return;

	n=get_descmax(zenodes,current_node);  // pour garder la couleur du nieme plus grand des groupes
	                                      // n est le numero du desc, pas son id

	color_clado( zenodes, zenodes[current_node].desc[n], whichcolor );


	for (i = 0; i < zenodes[current_node].nbdesc; i++)
		if(i!=n)
		{
			(*whichcolor)++;
			color_clado( zenodes, zenodes[current_node].desc[i], whichcolor );
	
		}
	
}

/*--------------------------------------------------*/
/*print recursively in memory the position of the nodes to be drawn
firts assign position to the leaves then going up assign position to the parents [ first_to_draw last_to_draw ]
*/
/*--------------------------------------------------*/
void print_clado(Node *zenodes, int current_node , FILE *fres, double echx, double echy, int sizestep, int s,int ss)
{
	int i;
	
	static int  feuilles = 0;

	if (zenodes[current_node].nbdesc == 0) //we have a leaf
	{
		zenodes[current_node].x = MARGECLADO; //x fixed
		zenodes[current_node].y = HAUTEURCOURBE + MARGECLADO + (feuilles * SIZEOFTEXT);
		zenodes[current_node].yy= feuilles * SIZEOFTEXT;
		zenodes[current_node].xx = 5;
		zenodes[current_node].first_to_draw=feuilles;

		feuilles++;
	}
	else
	{
		
		s = 0; //sum of all positions of desc
		ss=0;
		
		for (i = 0; i < zenodes[current_node].nbdesc; i++)//sea all the desc
		{
				
			zenodes[current_node].x = MARGECLADO + (zenodes[current_node].dist * echx);
			zenodes[current_node].xx= 5+ (zenodes[current_node].dist * echx);

			print_clado(zenodes, zenodes[current_node].desc[i], fres, echx, echy, sizestep, s,ss); // recursively call the fn
			s += zenodes[zenodes[current_node].desc[i]].y ;
			ss+= zenodes[zenodes[current_node].desc[i]].yy; 
			zenodes[current_node].first_to_draw=zenodes[zenodes[current_node].desc[i]].first_to_draw;


		}

		zenodes[current_node].y = (s / (float)zenodes[current_node].nbdesc);
		zenodes[current_node].yy = (ss / (float)zenodes[current_node].nbdesc);

	}


	return;

}

/*outputs a newick string from an array of struct Node*/
void nwkOut ( Node *zenodes, FILE *f ,int whichnode)
{
	int i;

	if (zenodes[whichnode].nbdesc==0)  {                                          

		fprintf(f,"%s_GR%d", zenodes[whichnode].name,zenodes[whichnode].specnumber);
	

		fprintf(f,":%lf",zenodes[whichnode].dist);

	}
	else {
		fprintf(f,"(");
		for (i=0;i<zenodes[whichnode].nbdesc;i++) {	
			nwkOut( zenodes, f ,zenodes[whichnode].desc[i]);
			if (i!=zenodes[whichnode].nbdesc -1)
				fprintf(f,",");
		}
		
			fprintf(f,")");
		
	
				fprintf(f,":%f",zenodes[whichnode].pval);
	}

}

/*--------------------------------------------------*/
/*Put a distance matrix (n  n) into a struct of size (n*(n-1))/2 which can be sorted according to distance*/

void mattolist(DistPair *dist_list , DistMat *d , float *max, float *min)
{
	int i,
	    j,
	    k = 0;

	*max = *min = d->dist[0][0];

	for (i = 0; i < d-> n - 1; i++)
		for (j = i + 1; j < d->n; j++)
		{
			dist_list[k].a = i;
			dist_list[k].b = j;
			dist_list[k].d = d->dist[i][j];
			
			if (*max < d->dist[i][j])
				*max = d->dist[i][j];
			if (*min > d->dist[i][j])
				*min = d->dist[i][j];

			dist_list[k].g = -1;

			k++;
		}

	qsort(dist_list, k, sizeof(DistPair), compareCase);


}





/*--------------------------------------------------*/
//END OF STATISTICS
/*--------------------------------------------------*/
/*************************************************/
int comparaison(const void *v1, const void *v2)
{
const float *fv1 = (float *)v1;
const float *fv2 = (float *)v2;
//printf("%f %f\n",*fv1,*fv2);
if (*fv1< *fv2)
	return(-1);
else 
	return (1);
}



/*--------------------------------------------------*/
//DRAWINGS
/*--------------------------------------------------*/




void createSVGhisto(char *ledir,struct DistanceMatrix dist_mat,int nbbids ,Results *scores,int nbresults,char *workdir)
{
int i,j,k;
int *histo;
float maxi=0,maxidist;
char chaine [12];

int largeur=360;
	int hauteur=260;

//	int largeur=720;
//	int hauteur=520;
	int marge=40;
	int bordure=60;
	int sizelegend=45;
	int x1,x2,y1,y2;

	float pas;
	FILE *svgout;

	double intervalle,echellex,echelley;
	char filename[256];
	float *histocum;
	char  *colors[3]={"#FFFFFF","#ffb83b","#00a5eb"};	
	int nbcomp=((dist_mat.n * (dist_mat.n -1)))/2.0;
	
		sprintf(filename,"%s%s/myHisto.svg",workdir,ledir); // all squares are pointing to a different file

	
	//sprintf(filename,"%s.disthist.svg",file);

	svgout=fopen(filename,"w");
	if (svgout==NULL)printf("no file in %s\n",filename),exit(1);
	fprintf(svgout,"<svg xmlns=\"http://www.w3.org/2000/svg\"  ");
	fprintf(svgout,"width=\"%d\" height=\"%d\" >\n",largeur+sizelegend, hauteur+sizelegend);

		
	histo=malloc(sizeof(int)*nbbids+1);
	if (histo==NULL)
	fprintf(stderr,"pb malloc histo(1)\n"),exit(1);
	
	histocum=malloc(sizeof(float)*nbcomp+1);
	if (histo==NULL)
	fprintf(stderr,"pb malloc histo(2)\n"),exit(1);
	
	for (i=0;i<nbbids;i++)histo[i]=0;

	k=0;
	for (i=0;i<dist_mat.n-1;i++)
		{
		for (j=i+1;j<dist_mat.n;j++)
			{
			if (maxi<dist_mat.dist[i][j])
				maxi=dist_mat.dist[i][j];
			histocum[k++] = (float) dist_mat.dist[i][j];
			
			}
	}	
	maxidist=maxi;
	qsort(histocum, nbcomp,sizeof(float),comparaison);
	

	intervalle=maxi/(float)nbbids;
	k=0;
	for (i=0;i<dist_mat.n-1;i++)
		for (j=i+1;j<dist_mat.n;j++)
			{
			k=dist_mat.dist[i][j]/intervalle;
			histo[k]++;
			}
			

	maxi=0;
	
for (i=0;i<nbbids;i++)
	{
	if (maxi<histo[i])
		maxi=histo[i];

	}
	fflush(stdout);
	largeur=largeur -bordure;
	hauteur=hauteur -bordure; 
		
	

	echellex=(float)largeur/nbbids;
	echelley=(float)hauteur/maxi;
	fprintf(svgout,"<g>\n");
	
	fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n",marge,marge,marge,hauteur+marge	 );
	fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n", marge ,hauteur+marge ,largeur+marge,hauteur+marge);

	pas=maxi/10.0;

	for(i=0;i<10;i++)
		{
			y1=hauteur+marge - ((i+1)*echelley*pas) ;
			fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n",marge-3, y1,marge,y1 );
			fprintf(svgout,"<text x=\"%d\" y=\"%d\" style=\"font-family: monospace; font-size: 8px;\">%d</text>\n",
			5,y1,(int)((i+1)*pas));	
			}
			
			
	fprintf(svgout,"<text x=\"5\" y=\"15\" style=\"font-family: monospace; font-size: 10px;\">Nbr</text>\n");
	  fprintf(svgout,"<text x=\"270\" y=\"280\" style=\"font-family: monospace; font-size: 10px;\">Dist. value</text>\n");

//plotting squares and values and ticks on x axis; write the image map using the exact same values  
	for (i=0;i<nbbids;i++)
		{
		//plot the value
		x1=marge+ ((i)*echellex);
		y1=hauteur+marge-(histo[i]*echelley);
		y2=(histo[i]*echelley) ;
		fprintf(svgout,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"  fill= \"%s\"  style=\" stroke: black;\" />\n", x1, y1, (int)echellex,y2,colors[1]);

		if ((nbbids<=20) || (nbbids>20 && i%2==0)) // because too much people on x axis if too much bids
	   		{
	   		sprintf(chaine,"%.2f",i*intervalle);
			fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n",x1,marge+hauteur, x1,marge+hauteur+5);
 			fprintf(svgout,"<text x=\"%d\" y=\"%d\" transform=\"rotate(90,%d,%d)\" style=\"font-family: monospace; font-size: 8px;\">%s</text>\n", 
					x1,marge+hauteur+6,x1,marge+hauteur+6,chaine);
	  		}
		
		}



echellex=largeur/maxidist;
for (j=0,i = nbresults - 1; i >= 0; i--)
		if ( scores[i].rank_general <=10) 
		{
//		x1=marge+(echellex*scores[i].d)	;
		x1=marge+(echellex*scores[i].d_jump)	;
		y1=hauteur+marge;
		y2=marge;
		
		fprintf(svgout,"<line id=\"CutDis%d\" visibility=\"hidden\" x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: red;\"/>\n",scores[i].rank_general,x1,y1, x1,y2);
	//	fprintf(svgout,"<line id=\"CutDis%d\" visibility=\"hidden\" x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: red;\"/>\n",j,x1,y1, x1,y2);

		j++;
		}


	fprintf(svgout,"</g>\n");
	fprintf(svgout,"</svg>\n");
	fclose(svgout);

	/**********************************AUTRE HISTO*/
	sprintf(filename,"%s%s/myHistoCum.svg",workdir,ledir); // all squares are pointing to a different file
	svgout=fopen(filename,"w");
	fprintf(svgout,"<svg xmlns=\"http://www.w3.org/2000/svg\"  ");
	fprintf(svgout,"width=\"%d\" height=\"%d\" >\n",largeur+sizelegend, hauteur+sizelegend+marge);
	fflush(stdout);
        maxi=histocum[nbcomp-1];
        echelley=(float)hauteur/maxi;   

       fprintf(svgout,"<g>\n");
        fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n",marge,marge,marge,hauteur+marge    );
        fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n", marge ,hauteur+marge ,largeur+marge,hauteur+marge);


        pas=hauteur/10;
        for(i=0;i<10;i++)
                {
                        
                        y1=(hauteur+marge)- ((i+1)*pas) ;
                        fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n",marge-3, y1,marge,y1);
                        sprintf(chaine,"%.2f",(float)(i+1)*(maxi/10));
                        k=8*(int)strlen(chaine);
                        x1= ((marge-k)>0) ?(marge-k):1;
                       fprintf(svgout,"<text x=\"%d\" y=\"%d\" style=\"font-family: monospace; font-size: 8px;\">%s</text>\n",x1, y1 ,chaine);
                       // fprintf(svgout,"<text x=\"%d\" y=\"%d\" transform=\"rotate(90,%d,%d)\" style=\"font-family: monospace; font-size: 8px;\">%s</text>\n", 
					//x1,marge+hauteur+5,x1,marge+hauteur+5,chaine);
                        }
                        
                        
        //drawing x axis
        echellex=(float)largeur/(float)nbcomp;
        for (i=0;i<10;i++)
                        {
                        k=(i+1)*((float)largeur/10.0);
                        x1=marge+ k;

                        sprintf(chaine,"%d",(int)((i+1)*(float)(nbcomp/10)));
                        fprintf(svgout,"<line x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: black;\"/>\n" , x1,marge+hauteur, x1,marge+hauteur+5);
                      	fprintf(svgout,"<text x=\"%d\" y=\"%d\" transform=\"rotate(90,%d,%d)\" style=\"font-family: monospace; font-size: 8px;\">%s</text>\n", 
                                x1,marge+hauteur+5,x1,marge+hauteur+5,chaine);
                        }
                
                fprintf(svgout,"<text x=\"%d\" y=\"%d\" style=\"font-family: monospace; font-size: 10px;\">Dist. value</text>\n",5,8);
	    	fprintf(svgout,"<text x=\"%d\" y=\"%d\" style=\"font-family: monospace; font-size: 10px;\">Rank</text>\n",320,230);
                        
        fprintf(svgout,"<polyline style=\"stroke: %s; stroke-width:1;fill: none;\"  points=\"",colors[2]); 
        x2=y2=0;
        for (i=0;i<nbcomp-1;i++)
                {
                        x1=marge+ ((i)*echellex);
                        y1=(hauteur+marge) -(histocum[i]*echelley) ;
                        if (i==0 || x1!=x2 || y1!=y2)
                                fprintf(svgout,"%d %d,",x1,y1); //draw new coords only
                        x2=x1;
                        y2=y1;
                        
                        }
        x1=marge+ ((i)*echellex);       
        y1=(hauteur+marge) -(histocum[i]*echelley)  ;
                fprintf(svgout,"%d %d\"/>",x1,y1);  
      //for   

	
	j=0;
for (i = nbresults - 1; i >= 0; i--)
	if ( scores[i].rank_general <=10) 
		{
		x1=marge;
		x2=largeur+marge;
//		y1=marge+hauteur-(scores[i].d * echelley);
		y1=marge+hauteur-(scores[i].d_jump * echelley);
		fprintf(svgout,"<line id=\"CutDis%d\" visibility=\"hidden\" x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: red;\"/>\n",scores[i].rank_general,x1,y1, x2,y1);
	
//		fprintf(svgout,"<line id=\"CutDis%d\" visibility=\"hidden\" x1=\"%d\" y1=\"%d\"  x2=\"%d\" y2=\"%d\" style=\" stroke: red;\"/>\n",j,x1,y1, x2,y1);
		j++;
		}



           
        
        fprintf(svgout,"\n</g>\n");
        fprintf(svgout,"</svg>\n");
        fclose(svgout);







	free(histo);


}



/*--------------------------------------------------*/

void CreateCurve2(Results *scores,  int nbresults, char *dirfiles, char *dataFilename, FILE *ff, double maxDist,  FILE *svgout,long nb_seq,double max_score,double min_score,int largeur,float minDi,float maxDi)
{
	int i,j=0, x = 0, y = 0, xo = 0, yo = 0,yo1,c,
	       bord = 20,
	       hauteur = HAUTEURCOURBE,
	      
	    
	       nbticks;

	float ytics,xtics;
	char *LavaColors[5]={"#320000","#F00000","#FFAF00","#FFFF00", "#CCCCCC"};
	double echellex,	echelleY1 ;
	char leg[1024];
	int hauteurMax=(nb_seq* SIZEOFTEXT)+HAUTEURCOURBE+(2*bord);
	if (ff == NULL)
		ff = stderr;
//axes
	

	echelleY1 = (double)hauteur / (double) (max_score);
//	echelleY2= (double)hauteur / (double) best_res[0][6];
	echellex= (double)largeur/(double) (scores[nbresults-2].d );



	svgFilledRectangleNoBorder(svgout, MARGECLADO+ (minDi*echellex), bord, ((maxDi-minDi)*echellex), hauteur,"lightblue");



	svgLine(svgout, MARGECLADO, bord, MARGECLADO, hauteur + bord, "black");
	svgString(svgout, 15, MARGECLADO, 15, "Asap Score", "black");
	//svgLine(svgout, MARGECLADO, hauteur + bord, largeur + MARGECLADO, hauteur + bord, "black");
		
//svgLine(svgout, MARGECLADO, hauteur + bord, (scores[nbresults-1].d *echellex)+ MARGECLADO, hauteur + bord, "black");
	
//svgString(svgout, 10, (scores[nbresults-1].d *echellex)+ MARGECLADO, hauteur + bord + 12, "dist", "black");

	nbticks = 10;

	ytics = max_score/ (float)nbticks;
	for (i = 0,j=10; i < nbticks; i++,j--)
	{
		y =  bord + (i*ytics*echelleY1) ;
		sprintf(leg, "%.1f", (i*ytics)+1);
		svgString(svgout, 15, MARGECLADO - 50, y, leg, "black");
		svgHorLine(svgout, MARGECLADO - 10, MARGECLADO, y, "black");
	}
	xtics=maxDist/10.0;
	for (i = 1; i <= 7; i++)
	{
		x= (i*xtics *echellex)+ MARGECLADO;
		sprintf(leg, "%.3f", (i*xtics));
		svgString(svgout, 15,x, hauteur+15+ bord, leg, "black");
		svgVertLine(svgout, x , hauteur+ bord,hauteur+10+ bord, "black");
	}

svgLine(svgout, MARGECLADO, hauteur + bord, x, hauteur + bord, "black");
	
svgString(svgout, 10,x-20, hauteur + bord + 12, "dist", "black");


ytics = nb_seq / (float)nbticks;

	xo = yo = yo1 =-1;

for (i = nbresults - 1,j=0; i >= 0; i--)

	{

		x = (double)(scores[i].d * echellex) + MARGECLADO;
		
		//y = hauteur + bord - ((scores[i].score) * echelleY1);
		y = hauteur + bord - ((max_score-scores[i].score) * echelleY1);
		if (xo != -1)
				{
		
					svgLine(svgout, xo, yo, x, y, "black");
		
				}
	if(scores[i].rank_general<=10)
			{
			sprintf(leg, "'inter/intra: %2.2e Nb species: %d / %d [%d groups] dist:%f'", scores[i].other_parameter, scores[i].nbspec, scores[i].nbspecRec, scores[i].nbgroups, scores[i].d);
	
			c=getCircleColor(scores[i].proba);
		//	fprintf(svgout, "    <line   visibility=\"hidden\"  id=\"bestline%d\"   x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"green\" />\n", j, x, HAUTEURCOURBE+bord, x, hauteurMax+bord);
		//fprintf(svgout, "    <line   visibility=\"hidden\"  id=\"bestline%d\"   x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"green\" />\n", scores[i].rank_general, x, HAUTEURCOURBE+bord, x, hauteurMax+bord);
			fprintf(svgout, "    <line   visibility=\"hidden\"  id=\"bestline%d\"   x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"green\" />\n", scores[i].rank_general, x, y, x, hauteurMax+(2*bord));
				
//fprintf(svgout, " <circle id=\"bestScore%d\" onmousemove=\"ShowTooltip(evt, %s)\"   onmouseout=\"HideTooltip(evt)\" onclick=\"circle_click(evt,%d,%f,%d,%d);\" cx=\"%d\" cy=\"%d\" r=\"7\" style=\"fill:%s;stroke:black;\" />\n", 
//				i,leg,j, scores[i].d, scores[i].nbspecRec,scores[i].rank_general, x, y,LavaColors[c]);
			fprintf(svgout, " <circle id=\"bestScore%d\" onmousemove=\"ShowTooltip(evt, %s)\"   onmouseout=\"HideTooltip(evt)\" onclick=\"circle_click(%d);\" cx=\"%d\" cy=\"%d\" r=\"7\" style=\"fill:%s;stroke:black;\" />\n", 
				i,leg,scores[i].rank_general, x, y,LavaColors[c]);
			j++;

			}

			
	
	
		xo = x;	
		yo = y;


	}

	draw_legend(svgout, 10, hauteur + 20);

}



/*--------------------------------------------------00a5eb;	 */
void write_javascript_svg(FILE *svgout)
{

	fprintf(svgout, "<script language=\"javascript\">\n");
	//fprintf(svgout, "var oldsvg;\n");
	fprintf(svgout, "var red_is_on=0\n");
	fprintf(svgout, "var oldy=-1;\n");
	

	fprintf(svgout, "function circle_click(n) {\n");
	fprintf(svgout, "		var svgDoc,svg_courbe,lalig;\n");
	fprintf(svgout, "      	a_effacer=0;\n");
	fprintf(svgout, "		if (oldy!=-1)\n");
	fprintf(svgout, "		{\n");
	fprintf(svgout, "		svgItemold = document.getElementById(\"bestline\"+oldy);\n");
	fprintf(svgout, "		svgItemold.setAttribute(\"visibility\",\"hidden\");\n\n");
	
	fprintf(svgout, "      	svg_courbe = document.getElementById(\"svg1\");\n");
	fprintf(svgout, "       svgDoc= svg_courbe.contentDocument;\n");
	fprintf(svgout, "	    lalig=svgDoc.getElementById(\"CutDis\"+oldy);\n");
	fprintf(svgout, "      	lalig.setAttribute(\"visibility\",\"hidden\");\n\n");
	fprintf(svgout, "      	svg_courbe = document.getElementById(\"svg2\");\n");
	fprintf(svgout, "       svgDoc= svg_courbe.contentDocument;\n");
	fprintf(svgout, "	    lalig=svgDoc.getElementById(\"CutDis\"+oldy);\n");
	fprintf(svgout, "      	lalig.setAttribute(\"visibility\",\"hidden\");\n");
 	fprintf(svgout, "      	if (oldy==n)\n");
 	fprintf(svgout, "      		a_effacer=1;\n");
 

	  fprintf(svgout, "		tablig=document.getElementById(\"row_res\"+oldy);\n");
	  fprintf(svgout, "		tablig.style.backgroundColor=\"transparent\";\n");


	fprintf(svgout, "		}\n");
	
	fprintf(svgout, "	svgItem = document.getElementById(\"bestline\"+n);\n");
	fprintf(svgout, "   if (a_effacer==0)\n");
	fprintf(svgout, "		svgItem.setAttribute(\"visibility\",\"visible\");\n");

	 fprintf(svgout, "	svg_courbe = document.getElementById(\"svg1\");\n");
	fprintf(svgout, "	svgDoc= svg_courbe.contentDocument;\n");
	fprintf(svgout, "	lalig=svgDoc.getElementById(\"CutDis\"+n);\n");
	fprintf(svgout, "   if (a_effacer==0)\n");
	fprintf(svgout, "		lalig.setAttribute(\"visibility\",\"visible\");\n");
	 fprintf(svgout, "	svg_courbe = document.getElementById(\"svg2\");\n");
	fprintf(svgout, "	svgDoc= svg_courbe.contentDocument;\n");
	fprintf(svgout, "	lalig=svgDoc.getElementById(\"CutDis\"+n);\n");
//	fprintf(svgout, "	lalig=svgDoc.getElementById(\"CutDis\"+n);\n");
	fprintf(svgout, "      	if (a_effacer==0)\n");
	fprintf(svgout, "	lalig.setAttribute(\"visibility\",\"visible\");\n");
	fprintf(svgout, "   tablig=document.getElementById(\"row_res\"+n);\n");
	fprintf(svgout, "    if (a_effacer==0)\n");
	fprintf(svgout, "   	tablig.style.backgroundColor=\"#ffb83b\";\n");


	fprintf(svgout, "	oldy=n;\n");
	fprintf(svgout, "}\n");

	fprintf(svgout, "\n");

	fprintf(svgout, " 	function init(evt)\n");
	fprintf(svgout, "	{\n");
	fprintf(svgout, "	    if ( window.svgDocument == null )\n");
	fprintf(svgout, "	    {\n");
	fprintf(svgout, "		svgDocument = evt.target.ownerDocument;\n");
	fprintf(svgout, "	    }\n");

	fprintf(svgout, "	    tooltip = svgDocument.getElementById('tooltip');\n");
	fprintf(svgout, "	    tooltip_bg = svgDocument.getElementById('tooltip_bg');\n");

	fprintf(svgout, "	}\n");
	
fprintf(svgout, "\n");

	fprintf(svgout, "	function ShowTooltip(evt, mouseovertext)\n");
	fprintf(svgout, "	{\n");
	
	fprintf(svgout, "		   var rect = svgDocument.getBoundingClientRect();\n");
   	fprintf(svgout, "	     var x = evt.clientX + svgDocument.scrollLeft - rect.left;\n");
    fprintf(svgout, "	    var y = evt.clientY + svgDocument.scrollTop - rect.top;\n");

	
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"x\",x+11);\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"y\",y+27);\n");
	fprintf(svgout, "	    tooltip.firstChild.data = mouseovertext;\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"visibility\",\"visible\");\n");

	fprintf(svgout, "	    length = tooltip.getComputedTextLength();\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"width\",length+8);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"x\",x+8);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"y\",y+14);\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"visibility\",\"visibile\");\n");
	fprintf(svgout, "	}\n");
	
fprintf(svgout, "\n");

	fprintf(svgout, "	function HideTooltip(evt)\n");
	fprintf(svgout, "	{\n");
	fprintf(svgout, "	    tooltip.setAttributeNS(null,\"visibility\",\"hidden\");\n");
	fprintf(svgout, "	    tooltip_bg.setAttributeNS(null,\"visibility\",\"hidden\");\n");
	fprintf(svgout, "	}\n");

	fprintf(svgout, "\n");
	
	fprintf(svgout, "	function select_species(evt,first,under)\n");
	fprintf(svgout, "		{\n");
	fprintf(svgout, "		var i=0;\n");
	
	fprintf(svgout, "		u=Number(under);\n");
	fprintf(svgout, "		if (red_is_on==0)\n");
	fprintf(svgout, "			{\n");
	
	fprintf(svgout, "			for (i=0;i<under;i++)\n");
	fprintf(svgout, "				{\n");
	
	fprintf(svgout, "				x=first-i;\n");
	fprintf(svgout, "				s=\"specie_\"+x;\n");	

	fprintf(svgout, "			    elt= svgDocument.getElementById(s);\n");
//		fprintf(svgout, "				alert (elt +\" \"+s);\n");
	fprintf(svgout, "				elt.setAttribute('fill', 'red');\n");

//	fprintf(svgout, "				i=i-1;\n");
	fprintf(svgout, "				}\n");
	fprintf(svgout, "			red_is_on=1;\n");
	fprintf(svgout, "			}\n");

	fprintf(svgout, "		else\n");
	fprintf(svgout, "			{\n");
	fprintf(svgout, "			for (i=0;i<under;i++)\n");
	fprintf(svgout, "				{\n");
	
	fprintf(svgout, "				x=first-i;\n");
	fprintf(svgout, "				s=\"specie_\"+x;\n");	

	fprintf(svgout, "			    elt= svgDocument.getElementById(s);\n");

	fprintf(svgout, "				elt.setAttribute('fill', 'black');\n");

//	fprintf(svgout, "				i=i+1;\n");
	fprintf(svgout, "				}\n");
		fprintf(svgout, "			red_is_on=0;\n");
	fprintf(svgout, "			}\n");

	fprintf(svgout, "		}\n");	
	
	fprintf(svgout, "function hideSVG(svgName) {\n");

  fprintf(svgout, "var style = document.getElementById(svgName).style.display;\n");
  fprintf(svgout, "if (style === \"none\")\n");
  fprintf(svgout, "  document.getElementById(svgName).style.display = \"block\";\n");
  fprintf(svgout, "else\n");
  fprintf(svgout, "  document.getElementById(svgName).style.display = \"none\";\n");
 fprintf(svgout, "}\n");

	
	
	
	fprintf(svgout, "</script>\n");
	
}
/*--------------------------------------------------*/
int getCircleColor(double v)
{

	int color;



	if (v < 0)
		color = 4;
 
	else if (v <= 0.001)
		color = 0;
	
	else if (v <= 0.05)	
		color=1;
	else if (v <= 0.1)
		color = 2;
		else
		color = 3; /*yellow*/

	return (color);
}
int getCircleColor2(double v)
{

	int color;



	if (v < 0)
		color = 4;
 
	else if (v >= 50)
		color = 0;
	
	else if (v >=10)	
		color=1;
	else if (v >=2)
		color = 2;
		else
		color = 3; /*yellow*/

	return (color);
}

/*--------------------------------------------------*/
void draw_legend(FILE *f, int x, int y)
{
	char *LavaColors[5]={"#320000","#F00000","#FFAF00","#FFFF00", "#CCCCCC"};


//	char *LavaColors[8] = /*foncé-->clair*/
//	{"#320000", "#780000", "#F00000", "#FF8700", "#FFAF00", "#FFD700", "#FFFF00" ,"#CCCCCC"};
	char *LavaLeg[5] = {"&lt;0.001", "&lt;0.05", "&lt;0.1", "&gt;0.1","N/A"};
	int i, xx = x;
	y+=50;
	svgString( f, 12, xx, y + 12, "Legend:", "black");
	xx+=70;

	for (i = 0; i < 5; i++)
	{
		svgCircle(f, xx, y, 6, LavaColors[i]);
		svgString( f, 12, xx-5, y + 15, LavaLeg[i], "black");
		xx = xx + 62;

	}
	svgHorLine(f,xx+10,xx+20,y + 15,"#ffb83b");
//	svgCircle(f, xx+20, y, 6, "green");
//	svgString( f, 12, xx+30, y + 15, "inter/intra", "black");
}
/*--------------------------------------------------*/
void draw_clado(Node *zenodes, FILE *fres, int nbnodes, int debut,int largeur)
{
	int i, k, j, ymax = 0;
	long int ymin = LONG_MAX;
	char leg[1224];
	int desc;
//	char *LavaColors[8] = /*foncé-->clair gris en dernier*/
//	{"#320000", "#780000", "#F00000", "#FF8700", "#FFAF00", "#FFD700", "#FFFF00", "#CCCCCC"};
//	double pi_intra,pi_inter;
	char *LavaColors[5]={"#320000","#F00000","#FFAF00","#FFFF00", "#CCCCCC"};
	
	
	
	for (i = debut; i <= nbnodes; i++)
	{

		ymin = LONG_MAX; ymax = 0;
		sprintf(leg, "%.2e", zenodes[i].pval);


		for (k = 0; k < zenodes[i].nbdesc; k++)
		{

			desc = zenodes[i].desc[k];

			svgHorLine(fres, zenodes[desc].x, zenodes[i].x, zenodes[desc].y , "black");

			if (ymin > zenodes[desc].y)
				ymin = zenodes[desc].y;
			if (ymax < zenodes[desc].y)
				ymax = zenodes[desc].y;
		}
		svgVertLine(fres, zenodes[i].x, ymin, ymax, "black");
	}

	svgHorLine(fres, zenodes[nbnodes].x, zenodes[nbnodes].x + ((largeur - 100) / zenodes[nbnodes].round), zenodes[nbnodes].y , "black");

	for (i = 0; i < debut; i++)
	{

		sprintf(leg, "%.12s", zenodes[i].name);

		
		fprintf(fres, "  <text x=\"%ld\" y=\"%ld\" id=\"specie_%d\" onmousemove=\"ShowTooltip(evt, '%s')\"   onmouseout=\"HideTooltip(evt)\" style=\"color:black;font-size :10;font-family = monospace;\" >%s</text>\n",
 						zenodes[i].x - MARGECLADO, zenodes[i].y, zenodes[i].first_to_draw, zenodes[i].name,leg);
			
		svgHorLine(fres, zenodes[i].x - 4, zenodes[i].x, zenodes[i].y, "black");
	}

	for (i = debut,j=0; i <= nbnodes; i++,j++)
	{
		int col = getCircleColor(zenodes[i].pval);
		if (zenodes[i].to_draw==0) col=4;	//tracer en gris les cerises
		if (zenodes[i].pval != 1)
			sprintf(leg, "'Proba: %2.2e dist:%f  nbgroups:%d node:%d nbunder:%d '", zenodes[i].pval, zenodes[i].dist, zenodes[i].nbgroups,i, zenodes[i].nb_under);
		else
			sprintf(leg, "'Proba: NA dist:%f  nbgroups:%d node:%d nb_under:%d'", zenodes[i].dist, zenodes[i].nbgroups,i, zenodes[i].nb_under);

fprintf(fres, "<rect id=\"Pval%d\"  onmousemove=\"ShowTooltip(evt, %s)\"   onmouseout=\"HideTooltip(evt)\"  \
			onmousedown=\"select_species(evt,%d,%d)\" x=\"%ld\" y=\"%ld\"  width=\"12\" height=\"12\"  style =\"fill:%s;stroke:black\"/> \n", 
	i,leg,zenodes[i].first_to_draw,zenodes[i].nb_under, zenodes[i].x-6, zenodes[i].y-6, LavaColors[col]);

//		fprintf(fres, " <circle id=\"Pval%d\"  onmousemove=\"ShowTooltip(evt, %s)\"   onmouseout=\"HideTooltip(evt)\"  onmousedown=\"select_species(evt,%d,%d)\" cx=\"%ld\" cy=\"%ld\" r=\"7\" fill=\"%s\" />\n", 
//		i,leg,zenodes[i].first_to_draw,zenodes[i].nb_under, zenodes[i].x, zenodes[i].y, LavaColors[col]);
	}
/*ADD what is a tooltip*/

	fprintf(fres, "<rect class=\"tooltip_bg\" id=\"tooltip_bg\"\n");
	fprintf(fres, "      x=\"0\" y=\"0\" rx=\"4\" ry=\"4\"\n");
	fprintf(fres, "      width=\"55\" height=\"17\" visibility=\"hidden\"/>\n");
 	fprintf(fres, " <text class=\"tooltip\" id=\"tooltip\"\n");
 	fprintf(fres, "     x=\"0\" y=\"0\" visibility=\"hidden\">Tooltip</text>\n");



}


/*--------------------------------------------------*/
void draw_square(Node *zenodes,int nodetodraw,int *spgraph,FILE *fgroups,int marge,int group,int x)
{
int y;
int heigth_todraw= (zenodes[nodetodraw].nb_under)*SIZEOFTEXT,icolor; 
int width_todraw=30;
//char *kokol[16]={"cornflowerblue","AQUA","TEAL","OLIVE","GREEN","LIME","YELLOW","ORANGE","RED","MAROON","FUCHSIA","PURPLE","PURPLE","LightSalmon","GRAY","SILVER"};

char *kokol[16]={"#666699","#0074D9","#7FDBFF","#39CCCC","#3D9970","#2ECC40","#01FF70","#FFDC00","#FF851B","#FF4136","#af1661","#F012BE","#B10DC9","#EFEFEF","#AAAAAA","#DDDDDD"};

//char *biColors[2]={"#ffb83b","#00a5eb"};//alternatively color groups with these 2 colors...

int spec_todrawlast= spgraph[zenodes[nodetodraw].first_to_draw] ;// donne la derniere espece 
			//	int colo=(spgraph[zenodes[nodetodraw].first_to_draw]+1)%16;
y=zenodes[spec_todrawlast].yy; // donne le coin sup gauche de la chaine de carractere de la derniere espece
y+=SIZEOFTEXT; // pour descendre dun cran sinon le dernier nest pas compris ds la boite
y-=heigth_todraw; //on remonte car la fn marche comme ca
icolor=	zenodes[nodetodraw].color %16;
			

svgFilledRectangle(fgroups, x, y+marge, width_todraw, heigth_todraw,kokol[icolor]);
//svgFilledRectangleWhiteBorder(fgroups, x, y+marge, width_todraw, heigth_todraw,"gray");
			
if (zenodes[nodetodraw].nb_under!=1)
	fprintf(fgroups, "  <text x=\"%d\" y=\"%d\">%d </text>",x+2,y+marge+SIZEOFTEXT,zenodes[nodetodraw].nb_under);
}

/*--------------------------------------------------*/
int nothing_good_under(Node *zenodes,int nodetodraw,double seuil)
{
int l;
for (l=0;l<zenodes[nodetodraw].nbdesc;l++)	
		{
			int n=	zenodes[nodetodraw].desc[l];
			if (zenodes[n].pval<seuil && zenodes[n].nb_under>1)	
				return (0);	
		}
return (1)	;	
}
/*--------------------------------------------------*/
void draw_Rec_square(Node *zenodes,int nodetodraw,int *spgraph,FILE *fgroups,int marge,int group,int x,double seuil)
{
int l;
if (zenodes[nodetodraw].nb_under==1 || (zenodes[nodetodraw].pval>seuil && nothing_good_under(zenodes,nodetodraw,seuil)))
 	draw_square(zenodes,nodetodraw,spgraph,fgroups,marge,group,x);	
else
	for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba

						draw_Rec_square(zenodes,zenodes[nodetodraw].desc[l],spgraph,fgroups,marge,0,x,seuil);		
}	
/*--------------------------------------------------*/
void draw_Rec_square_old(Node *zenodes,int nodetodraw,int *spgraph,FILE *fgroups,int marge,int group,int x,double seuil)
{
int l;
if (zenodes[nodetodraw].pval>seuil ||zenodes[nodetodraw].nb_under==1)
 	draw_square(zenodes,nodetodraw,spgraph,fgroups,marge,group,x);	
else
	for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
						draw_Rec_square(zenodes,zenodes[nodetodraw].desc[l],spgraph,fgroups,marge,0,x,seuil);		
}
/*--------------------------------------------------*/	
void tagg_rec(Node *zenodes,int nodetodraw,int *spec,double seuil)
{
	int i;

if (zenodes[nodetodraw].nbdesc==0)
	zenodes[nodetodraw].specnumber= *spec;
else
	for (i=0;i<zenodes[nodetodraw].nbdesc;i++)
		{
		
		tagg_rec(zenodes,zenodes[nodetodraw].desc[i],spec,seuil);
		if (zenodes[nodetodraw].pval<=seuil && zenodes[zenodes[nodetodraw].desc[i]].pval<=seuil) (*spec)++;
		}	

}
/*--------------------------------------------------*/	
		
void tagg_allSpec(Results *scores,Node *zenodes,int partition,double seuil)
{
int i,nbspec=0;;
int nodetodraw;
for (i=0;i<scores[partition].nbgroups;i++)
	{
	nodetodraw=scores[partition].listNodes[i];
	tagg_rec(zenodes,nodetodraw,&nbspec,seuil);
	nbspec++;
	}
}	
/*--------------------------------------------------*/	
void ecrit_esp_sous_node(Node *zenodes,int nodetodraw,FILE *f)
{
int i;
	if (zenodes[nodetodraw].nbdesc==0)
		fprintf(f,"%s ",zenodes[nodetodraw].name);	
	else
		for (i=0;i<zenodes[nodetodraw].nbdesc;i++ )
			ecrit_esp_sous_node( zenodes,zenodes[nodetodraw].desc[i], f);
}


void ecrit_rec_esp_sous_node(Node *zenodes,int nodetodraw,FILE *f,float seuil,int *j)
{
int l;
if (zenodes[nodetodraw].nb_under==1 ||  nothing_good_under(zenodes,nodetodraw,seuil))
 	{	fprintf(f,"\nGroup[ %d ] n: %d ;id: ",(*j)+1,zenodes[nodetodraw].nb_under);ecrit_esp_sous_node(zenodes,nodetodraw,f);	}
else
	for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba

						{ecrit_rec_esp_sous_node(zenodes,zenodes[nodetodraw].desc[l],f,seuil,j);	(*j)++;	}
}
/*--------------------------------------------------*/	
void ecrit_esp_cvs(Node *nodetodraw,int thenode,FILE *f,int sprec)
{
int i;
	if (nodetodraw[thenode].nbdesc==0)
		fprintf(f,"%s , %d \n",nodetodraw[thenode].name,sprec+1);	
	else
		for (i=0;i<nodetodraw[thenode].nbdesc;i++ )
			ecrit_esp_cvs(nodetodraw,nodetodraw[thenode].desc[i], f,sprec);
}

void ecrit_rec_esp_cvs(Node *zenodes,int nodetodraw,FILE *f,float seuil,int *j)
{
int l;
if (zenodes[nodetodraw].nb_under==1 ||  nothing_good_under(zenodes,nodetodraw,seuil))
 	{	ecrit_esp_cvs(zenodes,nodetodraw,f,*j);	}
else
	for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba

						{ecrit_rec_esp_cvs(zenodes,zenodes[nodetodraw].desc[l],f,seuil,j);	(*j)++;	}
}


/*void ecrit_fichier_texte_presque( char *dirfiles,int nbres,Node *zenodes,Results *scores,FILE *fres,float seuil)
{

	
	int i,j=0,k,l;
FILE *f;
FILE *f_l;
char nom_fic[256];

 
j=0;

for (k=0;k<nbres;k++)
		{
		sprintf(nom_fic,"%s/group_%d",dirfiles,k+1);	
			f=fopen (nom_fic,"w");
			sprintf(nom_fic,"%s/groupe_%d",dirfiles,k+1);	
			f_l=fopen (nom_fic,"w");
			if (f==NULL){fprintf(fres,"<H1>open text file %s pb\n</H1>",nom_fic);return;}
			fprintf(f,"Partition %d \nScore: %d\nProba: %e\n nb groups:%d (%d)\n",k+1,scores[k].rank_general,scores[k].proba,scores[k].nbspecRec,scores[k].nbgroups);
			fprintf(f,"------------------------------------------------------------\n");
			
		if (scores[k].rank_general<=10 )
			{
				j=0;
			for (i = 0; i <scores[k].nbgroups; i++)

				{
				int nodetodraw=scores[k].listNodes[i];
				fprintf(f,"\n\n--->(pval %e <?%e)  node%d under=%d\n",zenodes[nodetodraw].pval,seuil,nodetodraw,zenodes[nodetodraw].nb_under);
					
				if ( (zenodes[nodetodraw].pval==1 || zenodes[nodetodraw].pval>seuil || zenodes[nodetodraw].nb_under==1 )	)
					{
						fprintf(f,"\nGroup[ %d ] n: %d ;id: ",j+1,zenodes[nodetodraw].nb_under);
						ecrit_esp_sous_node(zenodes,nodetodraw,f);

						//ecrit_esp_sous_node(zenodes,nodetodraw,f_l);
						j++;
					}
						
				else	
				
					for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
						{
							int aa=	zenodes[nodetodraw].desc[l];
							fprintf(f,"\n*Group[ %d ] n: %d ;id: ",j+1,zenodes[aa].nb_under);
							ecrit_esp_sous_node(zenodes,aa,f,seuil);
							//ecrit_rec_esp_sous_node(zenodes,aa,f_l);
							j++;
						}
//
					
				}
		
	
		

			}	

			fclose(f);
			fclose(f_l);
		}
		
}
*/




void ecrit_fichier_texte( char *dirfiles,int nbres,Node *zenodes,Results *scores,FILE *fres,float seuil)
{

	
	int i,j=0,k,l;
FILE *f;
FILE *f_l;
char nom_fic[256];

 
j=0;

for (k=0;k<nbres;k++)
	
	

		{
		sprintf(nom_fic,"%s/group_%d",dirfiles,k+1);	
			f=fopen (nom_fic,"w");
			sprintf(nom_fic,"%s/groupe_%d",dirfiles,k+1);	
			f_l=fopen (nom_fic,"w");
			if (f==NULL){fprintf(fres,"<H1>open text file %s pb\n</H1>",nom_fic);return;}
			fprintf(f,"Partition %d \nScore: %d\nProba: %e\n nb groups:%d (%d)\n",k+1,scores[k].rank_general,scores[k].proba,scores[k].nbspecRec,scores[k].nbgroups);
			fprintf(f,"------------------------------------------------------------\n");
			
		if (scores[k].rank_general<=10 )
			{
				j=0;
			for (i = 0; i <scores[k].nbgroups; i++)

				{
				int nodetodraw=scores[k].listNodes[i];
				//fprintf(f,"\n\n--->(pval %e <?%e) k:%d node%d round: %d under=%d dist:%f %f/%f\n",
					//zenodes[nodetodraw].pval,seuil,k,nodetodraw,zenodes[nodetodraw].round,zenodes[nodetodraw].nb_under,zenodes[nodetodraw].dist,scores[k].d_jump,scores[k].d);
				

				if (zenodes[nodetodraw].pval==1 || zenodes[nodetodraw].pval>seuil || zenodes[nodetodraw].nb_under==1  )

					{
					if (zenodes[nodetodraw].dist==scores[k].d) //on est pile poil a lendroit ou il faut casser
					  {
						if ( nothing_good_under(zenodes,nodetodraw,seuil)) //on splite juste les groupes
							{
							for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
								{ 

									int aa=	zenodes[nodetodraw].desc[l];
									fprintf(f,"\nGroup[ %d ] n: %d ;id: ",j+1,zenodes[aa].nb_under);
									ecrit_esp_sous_node(zenodes,aa,f);
									ecrit_esp_cvs(zenodes,aa,f_l,j);
									j++;
								}
							}
						else
							{

							for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
								{
								int aa=	zenodes[nodetodraw].desc[l];
								//fprintf(f,"\n*Group[ %d ] n: %d ;id: ",j+1,zenodes[aa].nb_under);
								ecrit_rec_esp_sous_node(zenodes,aa,f,seuil,&j);
						
								ecrit_rec_esp_cvs(zenodes,aa,f_l,seuil,&j);
								j++;
								}

							}


						}//NOT PILEPOIL 
					else{
						fprintf(f,"\nGroup[ %d ] n: %d ;id: ",j+1,zenodes[nodetodraw].nb_under);
						ecrit_esp_sous_node(zenodes,nodetodraw,f);
						ecrit_esp_cvs(zenodes,nodetodraw,f_l,j);
						
						//ecrit_esp_sous_node(zenodes,nodetodraw,f_l);
						j++;
					}
					
					}// end of simple
						
				else	
				
					for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
						{
							int aa=	zenodes[nodetodraw].desc[l];
							//fprintf(f,"\n*Group[ %d ] n: %d ;id: ",j+1,zenodes[aa].nb_under);
							ecrit_rec_esp_sous_node(zenodes,aa,f,seuil,&j);
						
							ecrit_rec_esp_cvs(zenodes,nodetodraw,f_l,seuil,&j);
							j++;
						}
//
					
				} //end of i
		
			}//end of scorek.
		fclose(f);
			fclose(f_l);

			}	

			
		
		
}


/*--------------------------------------------------*/
/*Draw the groups the way Nico wants it...*/
void draw_nico(Node *zenodes, FILE *fgroups, int nbspecies,Results *scores,int nresult,double seuil,int nbest,int nbnodes,int largeur_clado)
{
	int i, k, j, l,maxLengthName=0,largeur_nico,hauteur;
	long x;
	int largeur_species, largeur_classes;
	int largeurlettre=7;
	char leg[1224];
	int desc;
	int marge=30;
	long ymin,ymax;
	int joli_y=SIZEOFTEXT/2;
	char *LavaColors[5]={"#320000","#F00000","#FFAF00","#FFFF00", "#CCCCCC"};
	int *spgraph;

 	spgraph=malloc(sizeof(int)*nbspecies);

	
	//first estimate lentgh max of specie
		for (i = 0; i < nbspecies; i++)
		{
			k=strlen(zenodes[i].name);
			if (k>maxLengthName)
				maxLengthName=k;
		}
	
	// evaluate size of graphic	
	//LArgeur of soecies name +	
	largeur_species=(maxLengthName*largeurlettre)+5;

	largeur_classes=	(60 *nbest);

	largeur_nico=largeur_species+largeur_clado+	largeur_classes+( 3*marge);
	hauteur=nbspecies*SIZEOFTEXT;

	svgImageCreate(fgroups,largeur_nico, 2*marge+(hauteur));

for (i = 0; i < nbspecies; i++)
	{

		sprintf(leg, "%s", zenodes[i].name);


		fprintf(fgroups, "  <text x=\"%ld\" y=\"%ld\"  id=\"specie_%d\"  style=\"color:black;font-size :15;font-family = monospace;\" >%s </text>\n",
 						zenodes[i].xx , 
 						zenodes[i].yy+marge+SIZEOFTEXT, 
 						
 						zenodes[i].first_to_draw, 

 						leg);
 	spgraph[zenodes[i].first_to_draw]=i;
 										


					//here node<todraw
	}

j=0;
x=0;
fprintf (fgroups,"<text x=\"10\" y=\"10\" style=\"color:black;font-size :12;font-family = monospace;\">Nb groups </text>\n");
fprintf (fgroups,"<text x=\"10\" y=\"25\" style=\"color:black;font-size :12;font-family = monospace;\">Score</text>\n");
for (k=0;k<nresult;k++)
	
		if (scores[k].rank_general<=10 )
		{
			//int col=0;
			x=	largeur_species+(50*j)+marge; 
			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :12px;font-family = monospace;\" >[%d] </text>\n",x ,SIZEOFTEXT/2,scores[k].nbspecRec);
			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :10px;font-family = monospace;\" >[%.1f] </text>\n",x ,SIZEOFTEXT+2,scores[k].score);
//			if(scores[k].proba <=seuil)
//			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :12;font-family = monospace;\" >[%.1e] </text>\n",x ,SIZEOFTEXT+2,scores[k].proba);
//		 	else
//		 		fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :10;font-family = monospace;\" >*%2.1f*</text>\n",x ,SIZEOFTEXT,100*scores[k].other_parameter);
			for (i = 0; i <scores[k].nbgroups; i++)
//		 	for (i = 0; i <scores[k].nbspec; i++)
				{
				int nodetodraw=scores[k].listNodes[i];

				if ((zenodes[nodetodraw].round !=k) &&(zenodes[nodetodraw].pval>seuil || zenodes[nodetodraw].nb_under==1 )	)
					draw_square(zenodes,nodetodraw,spgraph,fgroups,marge,i,x);	
				else	//LEQUEL pour le parametre truc... en fait split de celui sur leuql on est... 
					for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
						draw_Rec_square(zenodes,zenodes[nodetodraw].desc[l],spgraph,fgroups,marge,i*l,x,seuil);	
//
				
				}
		
			tagg_allSpec(scores,zenodes,k,seuil);//tagg all species with a group nbr corresponding to res
	
			j++;

		}	

	int arbreagauche=x+marge+50;



	for (i = nbspecies; i <= nbnodes; i++)
	{

		ymin = LONG_MAX; ymax = 0;

		for (k = 0; k < zenodes[i].nbdesc; k++)
		{

			desc = zenodes[i].desc[k];

		svgHorLine(fgroups, zenodes[desc].xx+arbreagauche, zenodes[i].xx+arbreagauche, zenodes[desc].yy+marge+joli_y , "black");

			if (ymin > zenodes[desc].yy)
				ymin = zenodes[desc].yy;
			if (ymax < zenodes[desc].yy)
				ymax = zenodes[desc].yy;
		}
		
		svgVertLine(fgroups, zenodes[i].xx+arbreagauche, ymin+marge+joli_y, ymax+marge+joli_y, "black");
	}

	svgHorLine(fgroups, zenodes[nbnodes].xx+arbreagauche,  zenodes[nbnodes].xx+arbreagauche+100, zenodes[nbnodes].yy+marge+joli_y , "black");

	for (i = 0; i < nbspecies; i++)
	{
	svgHorLine(fgroups, 5+arbreagauche - 4, 5+arbreagauche, zenodes[i].yy+marge+joli_y, "black");	
	}

	for (i = nbspecies; i <= nbnodes; i++)
	{
		int col = getCircleColor(zenodes[i].pval);
		if (zenodes[i].to_draw==0) col=4;	//tracer en gris les cerises		
		fprintf(fgroups, " <circle id=\"Pval%d\"  cx=\"%ld\" cy=\"%ld\" r=\"7\" fill=\"%s\" />\n", 
		i,zenodes[i].xx+arbreagauche, zenodes[i].yy+marge+joli_y, LavaColors[col]);
	}


fprintf(fgroups, "</svg>\n");
	
free(spgraph);
}



/*--------------------------------------------------*/
/*Draw the groups the way Nico wants it...*/
void draw_nico_old(Node *zenodes, FILE *fgroups, int nbspecies,Results *scores,int nresult,double seuil,int nbest,int nbnodes,int largeur_clado)
{
	int i, k, j, l,maxLengthName=0,largeur_nico,hauteur;
	long x;
	int largeur_species, largeur_classes;
	int largeurlettre=7;
	char leg[1224];
	int desc;
	int marge=30;
	long ymin,ymax;
	int joli_y=SIZEOFTEXT/2;
	char *LavaColors[5]={"#320000","#F00000","#FFAF00","#FFFF00", "#CCCCCC"};
	int *spgraph;

 	spgraph=malloc(sizeof(int)*nbspecies);

	
	//first estimate lentgh max of specie
		for (i = 0; i < nbspecies; i++)
		{
			k=strlen(zenodes[i].name);
			if (k>maxLengthName)
				maxLengthName=k;
		}
	
	// evaluate size of graphic	
	//LArgeur of soecies name +	
	largeur_species=(maxLengthName*largeurlettre)+5;
	largeur_classes=	(60 *nbest);

	largeur_nico=largeur_species+largeur_clado+	largeur_classes+( 3*marge);
	hauteur=nbspecies*SIZEOFTEXT;

	svgImageCreate(fgroups,largeur_nico, 2*marge+(hauteur));

for (i = 0; i < nbspecies; i++)
	{

		sprintf(leg, "%s", zenodes[i].name);


		fprintf(fgroups, "  <text x=\"%ld\" y=\"%ld\"  id=\"specie_%d\"  style=\"color:black;font-size :15;font-family = monospace;\" >%s </text>\n",
 						zenodes[i].xx , 
 						zenodes[i].yy+marge+SIZEOFTEXT, 
 						
 						zenodes[i].first_to_draw, 

 						leg);
 	spgraph[zenodes[i].first_to_draw]=i;
 										


					//here node<todraw
	}

j=0;
x=0;
fprintf (fgroups,"<text x=\"10\" y=\"10\" style=\"color:black;font-size :12;font-family = monospace;\">Nb groups </text>\n");
fprintf (fgroups,"<text x=\"10\" y=\"25\" style=\"color:black;font-size :12;font-family = monospace;\">Score</text>\n");
for (k=0;k<nresult;k++)
	
		if (scores[k].rank_general<=10 )
		{
			//int col=0;
			x=	largeur_species+(50*j)+marge; 
			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :12;font-family = monospace;\" >[%d] </text>\n",x ,SIZEOFTEXT/2,scores[k].nbspec);
			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :12px;font-family = monospace;\" >[%.1f] </text>\n",x ,SIZEOFTEXT+2,scores[k].score);
//			if(scores[k].proba <=seuil)
//			fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :12;font-family = monospace;\" >[%.1e] </text>\n",x ,SIZEOFTEXT+2,scores[k].proba);
//		 	else
//		 		fprintf(fgroups,"<text x=\"%ld\" y=\"%d\"   style=\"color:black;font-size :10;font-family = monospace;\" >*%2.1f*</text>\n",x ,SIZEOFTEXT,100*scores[k].other_parameter);
			for (i = 0; i <scores[k].nbgroups; i++)
//		 	for (i = 0; i <scores[k].nbspec; i++)
				{
				int nodetodraw=scores[k].listNodes[i];

				if ((zenodes[nodetodraw].round !=k) &&(zenodes[nodetodraw].pval>seuil || zenodes[nodetodraw].nb_under==1 )	)
					draw_square(zenodes,nodetodraw,spgraph,fgroups,marge,i,x);	
				else	//LEQUEL pour le parametre truc... en fait split de celui sur leuql on est... 
					for (l=0;l<zenodes[nodetodraw].nbdesc;l++) //Ici voir comment on peut dessiner ca recursivement si ya des noeuds dessous a grosse proba
						draw_Rec_square(zenodes,zenodes[nodetodraw].desc[l],spgraph,fgroups,marge,i*l,x,seuil);	
//
				
				}
		
			tagg_allSpec(scores,zenodes,k,seuil);//tagg all species with a group nbr corresponding to res
	
			j++;

		}	

	int arbreagauche=x+marge+50;



	for (i = nbspecies; i <= nbnodes; i++)
	{

		ymin = LONG_MAX; ymax = 0;

		for (k = 0; k < zenodes[i].nbdesc; k++)
		{

			desc = zenodes[i].desc[k];

		svgHorLine(fgroups, zenodes[desc].xx+arbreagauche, zenodes[i].xx+arbreagauche, zenodes[desc].yy+marge+joli_y , "black");

			if (ymin > zenodes[desc].yy)
				ymin = zenodes[desc].yy;
			if (ymax < zenodes[desc].yy)
				ymax = zenodes[desc].yy;
		}
		
		svgVertLine(fgroups, zenodes[i].xx+arbreagauche, ymin+marge+joli_y, ymax+marge+joli_y, "black");
	}

	svgHorLine(fgroups, zenodes[nbnodes].xx+arbreagauche,  zenodes[nbnodes].xx+arbreagauche+100, zenodes[nbnodes].yy+marge+joli_y , "black");

	for (i = 0; i < nbspecies; i++)
	{
	svgHorLine(fgroups, 5+arbreagauche - 4, 5+arbreagauche, zenodes[i].yy+marge+joli_y, "black");	
	}

	for (i = nbspecies; i <= nbnodes; i++)
	{
		int col = getCircleColor(zenodes[i].pval);
		if (zenodes[i].to_draw==0) col=4;	//tracer en gris les cerises		
		fprintf(fgroups, " <circle id=\"Pval%d\"  cx=\"%ld\" cy=\"%ld\" r=\"7\" fill=\"%s\" />\n", 
		i,zenodes[i].xx+arbreagauche, zenodes[i].yy+marge+joli_y, LavaColors[col]);
	}


fprintf(fgroups, "</svg>\n");
	
free(spgraph);
}
/*--------------------------------------------------*/
/*No use here but keep in case we want a tree some day
void print_TreeAsap(char *treestring, int nbesp, FILE *svgout, char *path)
{
	int height = nbesp * 20;s

	fprintf(svgout, "<html>\n");
	fprintf(svgout, "<head>\n");
	fprintf(svgout, "<script type=\"text/javascript\" src=\"%s/js/raphael-min.js\"></script> \n", path);
	fprintf(svgout, "<script type=\"text/javascript\" src=\"%s/js/jsphylosvg-min.js\"></script>\n", path);

	fprintf(svgout, "<script type=\"text/javascript\">\n");

	fprintf(svgout, "window.onload = function(){\n");

	fprintf(svgout, "var dataObject = { newick: '%s' };\n", treestring);

	fprintf(svgout, "phylocanvas = new Smits.PhyloCanvas(\n");

	fprintf(svgout, "dataObject,\n");

	fprintf(svgout, "'svgCanvas',\n");

	fprintf(svgout, "2000, %d\n", height);

	fprintf(svgout, " );\n");

	fprintf(svgout, " };\n");

	fprintf(svgout, " </script>\n");



	fprintf(svgout, "</head>\n");


	fprintf(svgout, "<body>\n");
	fprintf(svgout, "<div id=\"svgCanvas\"> </div>\n");

	fprintf(svgout, "</body>\n");
	fprintf(svgout, "</html>	\n");

}

*/