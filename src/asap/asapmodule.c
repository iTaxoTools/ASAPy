/*
  ASAPy - Assemble Species by Automatic Partitioning with ASAP
  Copyright (C) 2021  Patmanidis Stefanos

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
Extension Module for ASAP

Consider using multi-phase extension module initialization instead:
https://www.python.org/dev/peps/pep-0489/
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include "asap.h"
#include "oldfns.h"

#include "wrapio.h"

#define SIGN( a ) ( ( (a) > 0 )?1: ( ((a)==0)?0:-1)  )
static int Increase(const void *v1, const void *v2){  	return (int)SIGN( *((double *)v1) - *((double *)v2));  };
#undef SIGN

// Set var = dict[str], do nothing if key does not exist.
// On failure, sets error indicator and returns -1.
// Return 0 on success.
// Beware: parsing strings allocates memory!
int parseItem(PyObject *dict, const char *str, const char t, void *var) {

	PyObject *item;
	PyObject *value;
	switch(t){
		case 'b':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(int *)var = PyObject_IsTrue(item);
			if (*(int *)var == -1) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected boolean value for key '%s'", str);
				return -1;
			}
			break;
		case 'i':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(int *)var = (int) PyLong_AsLong(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected integer value for key '%s'", str);
				return -1;
			}
			break;
		case 'd':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(double *)var = (double) PyFloat_AsDouble(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected double value for key '%s'", str);
				return -1;
			}
			break;
		case 'f':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			*(float *)var = (float) PyFloat_AsDouble(item);
			if (PyErr_Occurred()) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected float value for key '%s'", str);
				return -1;
			}
			break;
		case 's':
			item = PyDict_GetItemString(dict, str);
			if (item == NULL) return 0;
			value = PyUnicode_AsEncodedString(item, "utf-8", "~E~");
			if (value == NULL) {
				PyErr_Format(PyExc_TypeError, "parseItem: Expected string value for key '%s'", str);
				return -1;
			}
			const char *bytes = PyBytes_AS_STRING(value);
			*(const char **)var = malloc(strlen (bytes) + 1);
			strcpy(*(const char **)var, bytes);
			Py_XDECREF(value);
			break;
		default:
			PyErr_Format(PyExc_TypeError, "parseItem: Unexpected type: %c", t);
			return -1;
		}
	return 0;
}

// Code imported from ABGD
char *Built_OutfileName( char *file ){

	char * bout;
	int ii;

	char *simplename;

	bout = ( strrchr(file,'/') == NULL )? file : strrchr(file,'/')+1;        /* either the begining or after the last '/' */

	ii = ( strchr(bout,'.')==NULL )? strlen(bout) : strchr(bout,'.')-bout ;  /* # of char before the first '.' */


	simplename=malloc(sizeof(char)*ii+1);

	strncpy(simplename,bout,ii);

	simplename[ii]='\0';

	return simplename;
}

static PyObject *
asap_main(PyObject *self, PyObject *args, PyObject *kwargs) {

	/* module specific */

	PyObject *dict = kwargs;
	PyObject *item;
	FILE *fres = NULL;


	/* a bunch of beautiful structures usefull for the code */

	DistMat mat;               /* The matrix distance */
	DistPair *ListDistance;    /* distance for all sequence pairs */
	Composante comp;           /* The whole graph that describes the species */
	Results *scores;
	Tabcompo *strucompo;       /* each elemnt store how many groups and how many sequences in each group */
	Node *zenodes;             /* Nodes of the hierarchical clusterting tree */
	Parameter asap_param;  		/*stuff for asap*/
	Spart *myspar;

	int i,
//		*grp,
//	    nb_pairs,
	    nbresults = 0,
	    firstpart,
//	    n_best=0,
	    color_ori=5,

	    seed_asap=-1,
	    *no_node;       // report for each sequence, its current node


	FILE *f_in,

	     *fgroups,
//	 	 *ffgroups,
	     *svgout;

	char *fout,
			 *dirfiles,
	     *dirfiles_default = ".",
			 *file_data,
			 *file_log,
	     *file_res,
//	      *nametree,
	     *fname,
	     *namegroups,
	     *simple_name;
	char *meth[5]={ "K80_Kimura","JC69_Jukes-Cantor","N93_Tamura-Nei" , "Simple_Dist"};


//	char nametree[512];

	time_t t1, t2, t3, t5;

	char c;
	const char *thedate = NULL;
	const char *thedate_default = "?";

	int imethode = 1, fmeg = 0, withallfiles = 0;//imethode1 for Jukes
	int last_node;

	float maxDist,
	      min,
	      ts_tv = 2.0;     /* default value for the trans/transv rates for Kimura 2-p */
	int nbBestAsap=10;
	double best_score, echx, echy,max_score,min_score;

	int widthKlado;
	//float seuil_pvalue=0.05;

	float minAsapDist=0.005,maxAsapDist=0.05;

	// struct stat     statbuf;
	// struct tm      *tm;
	struct stat st = {0};

	// stat(argv[0], &statbuf);
	// tm = localtime(&statbuf.st_mtime);
	// strftime(thedate,80,"%FT%T", tm); // 11/19/20 - 05:34PM

		/*
			init
		 */
		dirfiles = NULL;
		t1 = time(NULL);

		asap_param.pond_pente=0.1;
		asap_param.pond_score=0.5;
		asap_param.replicates=1000;
		asap_param.seuil_pvalue=0.001;
		//asap_param.ledir="";
		asap_param.fres=stderr;
		asap_param.lenSeq=600;
		asap_param.onlyspart=0;

	if (!PyArg_ParseTuple(args, "s", &file_data)) return NULL;

	f_in = fopen(file_data, "r");
	if (f_in==NULL) {
		PyErr_Format(PyExc_FileNotFoundError, "asap_main: Input file not found: '%s'", file_data);
		return NULL;
	}

	if (parseItem(dict, "out", 's', &dirfiles)) return NULL;
	if (!dirfiles) dirfiles = dirfiles_default;

	if (dirfiles[strlen(dirfiles)-1]!='/')
		strcat(dirfiles,"/");

	if (stat(dirfiles, &st) == -1) {
		mkdir(dirfiles, 0700);
	}

	/*
		Header
	*/
	fflush(stdout);
	fflush(stderr);
	fprintf(stderr, "/*\n");
	fprintf(stderr, "\tASAP (Agglomerate Specimens by Automatic Processing)\n");
	fprintf(stderr, "\twill delineate species in your dataset in a few moments.\n");
	fprintf(stderr, "\tRemember that the final cut remains yours!\n");
	fprintf(stderr, "*/\n");

	// Print these here so they are redirected if needed
	printf("> Passing parameters to ASAP\n", file_data);
	printf("- file = %s\n", file_data);
	printf("- dirfiles = %s\n", dirfiles);

	if (parseItem(dict, "time", 's', &thedate)) return NULL;
	if (!thedate) thedate = thedate_default;
	printf("- thedate = %s\n", thedate);

	if (parseItem(dict, "number", 'i', &nbBestAsap)) return NULL;
	printf("- nbBestAsap = %i\n", nbBestAsap);

	if (parseItem(dict, "method", 'i', &imethode)) return NULL;
	printf("- imethode = %i\n", imethode);

	if (parseItem(dict, "seed", 'i', &seed_asap)) return NULL;
	printf("- seed_asap = %i\n", seed_asap);

	if (parseItem(dict, "sequence_length", 'i', &(asap_param.lenSeq))) return NULL;
	printf("- asap_param.lenSeq = %i\n", asap_param.lenSeq);

	if (parseItem(dict, "replicates", 'i', &(asap_param.replicates))) return NULL;
	printf("- asap_param.replicates = %i\n", asap_param.replicates);

	if (parseItem(dict, "seuil_pvalue", 'f', &(asap_param.seuil_pvalue))) return NULL;
	printf("- seuil_pvalue = %f\n", asap_param.seuil_pvalue);

	if (parseItem(dict, "pond_pente", 'f', &(asap_param.pond_pente))) return NULL;
	printf("- pond_pente = %f\n", asap_param.pond_pente);

	if (parseItem(dict, "rate", 'f', &ts_tv)) return NULL;
	printf("- ts_tv = %f\n", ts_tv);

	if (parseItem(dict, "all", 'b', &withallfiles)) return NULL;
	printf("- withallfiles = %i\n", withallfiles);

	if (parseItem(dict, "mega", 'i', &fmeg)) return NULL;
	printf("- fmeg = %i\n", fmeg);


	printf("\n> Begin ASAP core:\n\n");
	fflush(stdout);
	fflush(stderr);

	if (seed_asap== -1)
		srand( time(NULL) );
	else
		srand(seed_asap);

	simple_name = Built_OutfileName(file_data);


	if (dirfiles == NULL)
	{
		dirfiles = (char *) malloc( (size_t) sizeof(char) * 3);
		if(!dirfiles)fprintf(stderr, "main: cannot allocate dirfiles bye\n"), exit(2);
		dirfiles[0] = '.'; dirfiles[1] = '/';dirfiles[2] = '\0';


	}

	/*

	*/
	fout = (char * )malloc( (size_t) sizeof(char) * (strlen(dirfiles) + strlen("partitions.tab") + 1));
	if (!fout)fprintf(stderr, "main: cannot allocate fout bye\n"), exit(2);

	namegroups=malloc(sizeof(char)*( (strlen (dirfiles) + 20)));
	sprintf(namegroups,"%sgroups.svg",dirfiles);
	sprintf(fout, "%spartitions.tab", dirfiles);

	asap_param.f_out = fopen(fout, "w+");
	asap_param.fres=stdout;
	asap_param.web=0;
	if (asap_param.f_out == NULL)fprintf(stderr,"cannot open the output file %s, bye\n", fout), exit(1);

	fname = (char *) malloc( (size_t) sizeof(char) * (strlen (dirfiles) + strlen (simple_name) +11) );
	sprintf(fname, "%sspecies.svg", dirfiles);

	svgout = fopen(fname, "w");
	if (svgout == NULL)fprintf(stderr, "cannot open the graphic output file %s, bye\n", fname), exit(1);


	/*
		Read or build the distance matrix
	*/
	c = fgetc(f_in);
	rewind(f_in);
	if ( c == '>'){
		fprintf(stderr, "> asap is reading the fasta file and computing the distance matrix\n");
		// mat = compute_dis(f_in, imethode, ts_tv, &(asap_param.lenSeq),"",stdout);
		mat = compute_dis(f_in, imethode, ts_tv, &(asap_param.lenSeq),asap_param);
	}
	else
	{
		fprintf(stderr, "assuming a distance matrix file. Reading the matrix\n");
		if (fmeg==0)
			mat = read_distmat(f_in, ts_tv, NULL, NULL);
		else
			if (fmeg==1)
			{
				c=fgetc(f_in);
				fputc(c,f_in);

				if (c==',')
				{
					read_mega10(f_in,&mat);printf("done 10\n");
				}
				else
				readMatrixMegaCVS(f_in,&mat);
			}
		//else
			//readMatrixMega(f_in,&mat);
	}
	fclose(	f_in);
	fprintf(stderr,"End of matrix distance\n");
	myspar=malloc(sizeof(Spart)*mat.n);

		/*for (i=0;i<mat.n;i++)
			{
				//myspar[i].name=malloc(sizeof(char)*strlen( mat.names[i])+1);
				//strcpy(myspar[i].name,mat.names[i]);

				myspar[i].specie=malloc(sizeof(int)* nbBestAsap);

			}*/

	if (mat.n<MAXSPECIESGRAPH)
	widthKlado=WIDTHCLADO/3;
	else
	widthKlado=WIDTHCLADO;
	fprintf(stderr, "  %ld input sequences\n", mat.n);

	t2 = time(NULL);

	/*
		Get memory for needed struct
	*/

	asap_param.nbpairs = (mat.n * (mat.n - 1)) / 2;

	ListDistance = (DistPair *) malloc( (size_t) sizeof(DistPair) *  asap_param.nbpairs);
	if (!ListDistance)fprintf(stderr, "main: cannot allocate  ListDistance bye\n"), exit(2);
//
	no_node = (int *)malloc( (size_t) sizeof(int) * mat.n);               //indice qui me donne pour une seuqnec i a quel noeud elle est liée à chaque etape de lagglutnage
	if (!no_node)fprintf(stderr, "main: MEMORY ERROR error can allocate nonode bye\n"), exit(2);

	zenodes = (Node *) malloc( (size_t) sizeof(Node) * ((mat.n*2)-1));
	if (!zenodes)fprintf(stderr, "main: MEMORY ERROR error can allocate  zenodes bye\n"), exit(2);

	strucompo = (Tabcompo *) malloc( (size_t) sizeof(Tabcompo) * mat.n);
	if (!strucompo)fprintf(stderr, "main: cannot allocate  strucompo bye\n"), exit(2);


	scores = (Results *) malloc(  (size_t) sizeof(Results) * mat.n); //not enough juste a first dim reajusted when nbresulkts greater than mat.n
	if (!scores)fprintf(stderr, "main: cannot allocate  scores bye\n"), exit(2);
	//build the sorted struct from smallest dist to higher dist
	for (i=0;i<mat.n;i++)
				scores[i].listNodes=malloc(sizeof(int)*mat.n);

	initcomp(&comp, mat.n, stderr, "");
	inittabcompo(strucompo, mat.n, stderr, "");       /* the structures are oversized currently */
	initNodes(stderr, zenodes, mat, "");

	/*
		Set the first n nodes to their id --the leaves--
	*/
	for (i = 0; i < mat.n; i++)
		no_node[i] = i;


	/*
		from the distance matrix, build a sorted list of pairwise_distance, min and max
	*/
	mattolist(ListDistance , &mat , &maxDist, &min);


	//	for (i=0;i<nb_pairs;i++)
	//		fprintf(stderr,"%d %f %d %d\n",i,ListDistance[i].d,ListDistance[i].a,ListDistance[i].b);


	nbresults = 0;

	last_node = mat.n - 1;

	/*
		Run ASAP core
	*/

	fprintf(stderr,"> asap is building and testing all partitions\n  ");

// int do_agglutine(DistMat mat, Composante *comp, DistPair *ListDist, Results *scores, Tabcompo *strucompo, int nb_pairs, FILE *f_out,double *best, int *fi, FILE *ff, Node *zenodes, int *list_node, int *lastnode, char *ledir, int lenSeq, int replicates,float seuil_pvalue,float pond_pente)


	//nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo, nb_pairs, f_out, &best_score, &firstpart, stderr, zenodes, no_node, &last_node, "", len_seq, replicates,seuil_pvalue,pond_pente);
	nbresults = do_agglutine( mat, &comp, ListDistance, scores, strucompo,  &best_score, &firstpart,  zenodes, no_node, &last_node,asap_param);

fprintf(stderr,"> asap has finished building and testing all partitions\n  ");
		/*if (fdeb!=NULL)
		{
		fprintf(fdeb,"%d res\n",nbresults);
		for (i=0;i<nbresults;i++)
			fprintf(fdeb,"%d %d %d\n",i,scores[i].nbspec,scores[i].nbspecRec);
		fclose (fdeb);
	}*/

	qsort(scores,nbresults,sizeof (Results ),compareProba);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_proba=i+1;

	qsort(scores,nbresults,sizeof (Results ),compareParameter);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_pente=i+1;


	max_score=0;
	scores[0].score=(scores[0].rank_pente*(1.0-asap_param.pond_score))+(scores[0].rank_proba*asap_param.pond_score);
	min_score=scores[0].score;
	for (i = 0; i < nbresults+1; i++)
	{
		scores[i].score=(scores[i].rank_pente*(1.0-asap_param.pond_score))+(scores[i].rank_proba*asap_param.pond_score);
		if (max_score <scores[i].score)
			max_score =scores[i].score;
		if (min_score >scores[i].score)
			min_score =scores[i].score;
	}
	qsort(scores,nbresults,sizeof (Results ),compareRang);
	for (i = 0; i < nbresults+1; i++)
		scores[i].rank_general=i+1;



	file_res = (char * )malloc( (size_t) sizeof(char) * (strlen(dirfiles) + strlen("scores.tab") + 1));
	if (!file_res)fprintf(stderr, "main: cannot allocate file_res bye\n"), exit(2);
        sprintf(file_res,"%sscores.tab",dirfiles);
	fres = fopen(file_res,"w");

	fprintf(stderr, "> writing scores file\n",asap_param.lenSeq);
	fprintf(fres, "\n> %d Best asap scores (probabilities evaluated with seq length:%d)\n",nbBestAsap,asap_param.lenSeq);
	fprintf(fres, "  distance  #species   #spec w/rec  p-value pente asap-score\n");
	int nb_B=(nbresults<nbBestAsap)?nbresults:nbBestAsap;
	for (i = 0; i < nb_B; i++)

			{
						char toStar=' ';
					if (scores[i].d_jump>=minAsapDist && scores[i].d_jump<=maxAsapDist)
						toStar='*';

			 fprintf(fres, "%c%8.4f %8d  %12d  %.3e %e \t%f \n",
			    	toStar,
			 		scores[i].d_jump,
			 		scores[i].nbspec,
			      	scores[i].nbspecRec,
			       	scores[i].proba,

			       	 scores[i].other_parameter *100,
			       	 scores[i].score);
			}
	fclose(fres);

	/*if (withallfiles)
		//ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue);
	ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue,myspar,mat.n);*/
	//printf("creating histo in %s\n",dirfiles);
	//createSVGhisto(dirfiles,mat,20,scores, nbresults,WORKDIR_CL);
	if (asap_param.onlyspart==0)
		createSVGhisto(dirfiles,mat,20,scores, nbresults,"",simple_name);

	/*
		That
	*/
	if (asap_param.onlyspart==0)
	fprintf(stderr, "> asap is creating text and graphical output\n");
	else
		fprintf(stderr, "> asap is creating spart output\n");
	qsort(scores,nbresults,sizeof (Results ),compareSpecies);

if (asap_param.onlyspart==0)
	{
	fprintf(svgout, "<svg xmlns=\"http://www.w3.org/2000/svg\" onload=\"init(evt)\" ");
	fprintf(svgout, "width=\"%d\" height=\"%ld\" >\n", widthKlado + MARGECLADO + 20, HAUTEURCOURBE + MARGECLADO + ( mat.n * SIZEOFTEXT));


	CreateCurve2(scores, nbresults, dirfiles, simple_name, NULL, maxDist,  svgout,mat.n,max_score,min_score,widthKlado,minAsapDist,maxAsapDist);
	}
	clearalltab(strucompo, &comp, mat.n);



	echy = mat.n * SIZEOFTEXT;
	echx = widthKlado / (float)maxDist;

	print_clado(zenodes, last_node, NULL, echx, echy, (widthKlado - 100) / zenodes[last_node].round, 0,0);

	draw_clado(zenodes, svgout, last_node, mat.n,widthKlado);

	color_clado(zenodes, last_node,&color_ori);


//	draw_bestlines(zenodes, last_node,svgout, scores, echx, MARGECLADO + ( mat.n * SIZEOFTEXT) + HAUTEURCOURBE, nbresults,min_pvalue);

	fprintf(svgout, "</svg>\n");
	fclose(svgout);
	fgroups=fopen(namegroups,"w");
	if (fgroups==NULL)
	printf("Cant write %s\n",namegroups);
	else
	draw_nico(zenodes, fgroups, mat.n,scores,nbresults,asap_param.seuil_pvalue,10,last_node,widthKlado);
	//printf("write \n");
		for (i=0;i<mat.n;i++)
			{
				int n=zenodes[i].first_to_draw; //assign ed in print_clado
				if (n>=mat.n)printf("error***** %d \n",n);
				myspar[n].name=malloc(sizeof(char)*strlen( mat.names[i])+1);
				strcpy(myspar[n].name,mat.names[i]);
				//printf("i:%d %d %s %s\n",i,n,myspar[n].name,mat.names[i]);
				myspar[i].specie=malloc(sizeof(int)* nb_B+1);
				//myspar[i].specie_ori=malloc(sizeof(int)* nb_B+1);
				//myspar[i].group=malloc(sizeof(int)*(n+1));
			}
			/*for (i=0;i<mat.n;i++)
			printf("i:%d %s\n",i,myspar[i].name);*/
		//qsort(scores,nbresults,sizeof (Results ),compareRang);
	if (withallfiles)
	{
		int **o_sp;
			//ecrit_fichier_texte( dirfiles,nb_B, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue);
		qsort(scores,nbresults,sizeof (Results ),compareRang);


		o_sp=malloc(sizeof(int*)*nb_B);
		for (i=0;i<nb_B;i++)
			o_sp[i]=malloc(sizeof(int)*2);
		//fprintf(stderr,"go ecrit\n");
		ecrit_fichier_texte( dirfiles,nb_B-1,nbresults, zenodes,scores,asap_param.fres,asap_param.seuil_pvalue,myspar,mat.n,last_node,simple_name,asap_param.onlyspart);
		//fprintf(stderr,"go order\n");
		order_spart(o_sp,nb_B,myspar,mat.n);
		//fprintf(stderr,"go create\n");
		CreateSpartFile(myspar,dirfiles,nb_B,simple_name,stdout,scores,mat.n,thedate,"",meth[imethode],o_sp);
		//fprintf(stderr,"go print_Spart\n");
		//print_spart(myspar,nb_B,mat.n);
		CreateXMLFile(myspar,dirfiles,nb_B,simple_name,stdout,scores,mat.n,thedate,"",meth[imethode],o_sp);
		for (i=0;i<nb_B;i++)
			free(o_sp[i]);
		free(o_sp);
	}

	fclose(asap_param.f_out);
	fclose(fgroups);

	fprintf(stderr, "> results were write\n");

	t5 = time(NULL);

	resetcomp(&comp, mat.n);
	// fprintf(stderr, "  partition results are logged in: \n");
	// fprintf(stderr, "\tFull results: %s\n",fout);
	// fprintf(stderr, "\tThe csv file of rank x: %sgroupe_x\n",dirfiles);
	// fprintf(stderr, "\tThe result file of rank x: %sgroup_x\n",dirfiles);
	// fprintf(stderr, "\tThe graphic outputs are: %s*.svg\n",dirfiles);
	// fprintf(stderr, "\tSpart file is: %s%s.spart\n", dirfiles,simple_name);

	free (fout);
//	free(fileNex);
//	free(newickStringOriginal);
//	free(newickString);
	free (ListDistance);
	freecomp(&comp, mat.n);
	free_distmat(mat);



	t3 = time(NULL);

	fprintf(stderr, "> asap computation times were:\n");

	fprintf(stderr,"  %2ldm %2lds to read file and compute distance\n", (t2 - t1) / 60, (t2 - t1) % 60);
	fprintf(stderr,"  %2ldm %2lds to compute and test all partitions\n", (t3 - t2) / 60, (t3 - t2) % 60);
//	fprintf(stderr,"  %2ldm %2lds to build an nj tree\n", (t5 - t4) / 60, (t5 - t4) % 60);
	fprintf(stderr,"  --------------\n");
	fprintf(stderr,"  %2ldm %2lds total\n", (t5 - t1) / 60, (t5 - t1) % 60);

//	fclose(	f_out);

	Py_INCREF(Py_None);
	return Py_None;
}


static PyMethodDef AsapMethods[] = {
  {"main",  asap_main, METH_VARARGS | METH_KEYWORDS,
   "Run ASAP on a file for given parameters."},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyDoc_STRVAR(asap_doc,
"Assemble species by automatic partitioning.");

static struct PyModuleDef asapmodule = {
  PyModuleDef_HEAD_INIT,
  "asap",   /* name of module */
  asap_doc, /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module,
               or -1 if the module keeps state in global variables. */
  AsapMethods
};

PyMODINIT_FUNC
PyInit_asap(void)
{
	PyObject *m = NULL;

	if (!(m = PyModule_Create(&asapmodule)))
		return NULL;

	if (wrapio_init(m)) {
		Py_XDECREF(m);
		return NULL;
	}

	return m;
}
