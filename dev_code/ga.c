//Genetic algorithm code
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define num 1372        //number of atoms
#define tot 20
#define width 33        //width-1 is the width of coord file
#define width_ele 2     //width associated with atom type
#define tot_gen 3
#define rb 8.314e-3
#define T 300
//Roulette wheel parameters
#define div_fact 2 
#define pop_div 5
//checkpoint parameters
#define check 1

int randominlimits(int,int,int); //function to generate random number between two numbers.
void parent_gen(int []); //function to generate the parent reserviour.
int swap(char [][width],char [][width],char [][width_ele],char [][width_ele],int,double [],double []); //function to swap between two configurations and generate 2 CONFIG files.
void status_write(int);
int countlines(char *);
void ascend(double [],int);
int degen_count(double [],int);
void par_sel(int,double [],int []);
int chooseornot(int,double[]);
int child_accept(double);
void checkpoint(int,char *[],char *[]);

int main(int argc, char** argv)
{
	char *restart1; //first argument, containing retart instructions
	char *restart2; //second argument, containing the generation number
	int restart_val = -1;
	int restart_gen = -1;
	restart1 = *(argv+1);
	restart2 = *(argv+2);
	restart_val = atoi(restart1); //conversion of restart (char *) to int
	restart_gen = atoi(restart2);
	printf("%s\t%s\t%d\t%d\n",restart1,restart2,restart_val,restart_gen);
        if((restart_val != 1) && (restart_val != 0))
	{	
		printf("Enter 0 if you want fresh start and 1, if if it is restarting\n");
		exit(1);
	}
	else if((restart_val != 0) && (restart_gen == 1))
	{
		printf("Restarting the calculation and yet generation value not provided!\n");
		exit(1);
	}
	//cleaning the directory, to avoid appending on earlier files.
        //system("./clean.bash");
	srand(clock());

	int i,j,k,p,counter,u,v;
	int total = num*tot;
	int ret_const = -1;	// this variable would return a constant, whose value can be 1,2,3 and 4 depicting child 1, child 2, both children and none of children, respectively has/have been chosen. 
	FILE *f_atom;
	FILE *f_coord;
	char atfile[20]; //string to store the name of atomic name's file.
        char cordfile[20]; //string to store the name of coordinates of atom's file.
	int res[tot]; //array to store the identity of parents chosen.
	int gen = 0;
	char pos1[num][width];
        char pos2[num][width];
        char atom1[num][width_ele];
        char atom2[num][width_ele];
	double en[tot]; //array to store the energy of children, which are accepted to be part of population.
	double en_par[tot]; //array to store the energy of the parents.
	double acc_en[2]; //this array will store the energy of two children, which have been generated in particular GA cycle and have been accepted. Note that in each cycle there are 100 accpted children.
	double par_en[2]; //this array will store the energy of two parents, which will be swapped to generate two children.
	//dynamically allocating memory to store values of coordinates and atoms
	char *pos[total];
	for(i=0;i<total;i++)
	{
		pos[i] = (char *)malloc(width * sizeof(char));
		if(pos[i] == NULL)
		{
			printf("Issue with array allocation- 1\n");
			exit(1);
		}
	}	
	char *atom[total];
        for(i=0;i<total;i++)
        {
                atom[i] = (char *)malloc(width_ele * sizeof(char));
		if(atom[i] == NULL)
                {
                        printf("Issue with array allocation- 2\n");
                        exit(1);
                }
        }
	char *pos_inp[total];
	for(i=0;i<total;i++)
	{
		pos_inp[i] = (char *)malloc(width * sizeof(char));
		if(pos_inp[i] == NULL)
                {
                        printf("Issue with array allocation- 3\n");
                        exit(1);
                }
	}
	char *atom_inp[total];
	for(i=0;i<total;i++)
	{
		atom_inp[i] = (char *)malloc(width_ele * sizeof(char));
		if(atom_inp[i] == NULL)
                {
                        printf("Issue with array allocation- 4\n");
                        exit(1);
                }
	}

	char *pos_pre[total];
        for(i=0;i<total;i++)
        {
                pos_pre[i] = (char *)malloc(width * sizeof(char));
                if(pos_pre[i] == NULL)
                {
                        printf("Issue with array allocation- 5\n");
                        exit(1);
                }
        }
        char *atom_pre[total];
        for(i=0;i<total;i++)
        {
                atom_pre[i] = (char *)malloc(width_ele * sizeof(char));
                if(atom_pre[i] == NULL)
                {
                        printf("Issue with array allocation- 6\n");
                        exit(1);
                }
        }
	//reading files.
	counter = 0;
	for(p=0;p<tot;p++)
	{
		if(restart_val == 0)
			sprintf(atfile,"atomlist.%d",p);
		else if(restart_val == 1) //if, calculation is being restarted
			sprintf(atfile,"alist.backup.%d",p);
                else
			;	
        
        	f_atom = fopen(atfile,"r");
        	if(f_atom == NULL)
        	{
                	printf("Error in opening file-atomlist\n");
                	exit(1);
        	}
        	
                for(i = (0 + counter);i < (num + counter);i++)
        	{
                	for(j=0;j<width_ele;j++)
                	{
                        	fscanf(f_atom,"%c",&atom[i][j]);
                	}
        	}
        	fclose(f_atom);
        
        	if(restart_val == 0)
			sprintf(cordfile,"coord.%d",p);
		else if(restart_val == 1) //if, calculation is being restarted.
			sprintf(cordfile,"cord.backup.%d",p);
		else
			;
                f_coord = fopen(cordfile,"r");
        	if(f_coord == NULL)
        	{
                	printf("Error in opening file-coord\n");
                	exit(1);
        	}
        
        	for(i = (0 + counter);i < (num + counter);i++)
        	{
                	for(j=0;j<width;j++)
                	{
                        	fscanf(f_coord,"%c",&pos[i][j]);
                	}
        	}
        	fclose(f_coord);	
		//increasing counter
		counter = (counter + num);
        		
	}
	//choosing mating partners or parents.
	parent_gen(res);        
	//Generating input files for swap
	for(i=0;i<tot;i++)
	{
		u = *(res+i);
		for(j=0;j<num;j++)
		{
			for(k=0;k<width;k++)
			{
				pos_inp[((i * num) + j)][k] = pos[((num * u) + j)][k];
			}
		}
	}
	for(i=0;i<tot;i++)
        {
                u = res[i];
                for(j=0;j<num;j++)
                {
                        for(k=0;k<width_ele;k++)
                        {
                                atom_inp[((i * num) + j)][k] = atom[((num * u) + j)][k];
                        }
                }
        }
	
	//en_par[] array should be populated here.
	for(i=0;i<tot;i++) 
        {
		en_par[i] = 0.0; //may be remove this
        }
	//GA cycle starts here.
	//restart_gen = 1 for new calculations.
        for(gen=restart_gen;gen<=tot_gen;gen++) // TODO change to tot_gen here.
	{
		//swapping between two parents to generate two children.
		for(i=0;i<tot;i++) //initialising en[] array, which will store the energy of accepted children.
		{
			en[i] = 0.0;
		}
		k = 0;
		do 
		{
		 	acc_en[0] = 0.0;
			acc_en[1] = 0.0;
			par_en[0] = 0.0;
			par_en[1] = 0.0;
			//In case mating or swapping takes place between 2 randomly chosen mating partners from parent list. This might be relevant when acceptance criteria for children selection would be introduced.
  			u = rand() % tot;
			v = randominlimits(0,tot,u);
			
			par_en[0] = en_par[u];
			par_en[1] = en_par[v];
			 
			for(i=0;i<num;i++)
			{
				for(j=0;j<width;j++)
				{
					pos1[i][j] = pos_inp[((u * num) + i)][j];
					pos2[i][j] = pos_inp[((v * num) + i)][j];
				}
			}
			for(i=0;i<num;i++)
        		{
                		for(j=0;j<width_ele;j++)
                		{
                        		atom1[i][j] = atom_inp[((u * num) + i)][j];
                        		atom2[i][j] = atom_inp[((v * num) + i)][j];
                		}
        		}
	
			ret_const = swap(pos1,pos2,atom1,atom2,gen,acc_en,par_en);
			//checking, remove this
			FILE *fc;
			fc = fopen("check_stat.csv","a");
			fprintf(fc,"%d\t%d\t%lf\t%lf\t%d\n",gen,ret_const,acc_en[0],acc_en[1],k);
			fclose(fc);
			//This section stores only those configurations, which have been chosen/selected. If thr selection criteria of children has not been employed, both the configurations are stored for next generation.
			if(ret_const == 1) //if only child-1 is accepted.
			{
				for(i=0;i<num;i++) 
                                {
                                        for(j=0;j<width;j++)
                                        {
                                                pos_pre[((k * num) + i)][j] = pos1[i][j];
                                        }
                                        for(j=0;j<width_ele;j++)
                                        {
                                                atom_pre[((k * num) + i)][j] = atom1[i][j];
                                        }
                                }
				en[k] = acc_en[0];
				k = k + 1;
			}
			else if(ret_const == 2) //if only child-2 is accepted.
			{
				for(i=0;i<num;i++)
                                {
                                        for(j=0;j<width;j++)
                                        {
                                                pos_pre[(((k+1) * num) + i)][j] = pos2[i][j];
                                        }
                                        for(j=0;j<width_ele;j++)
                                        {
                                                atom_pre[(((k+1) * num) + i)][j] = atom2[i][j];
                                        }
                                }
				en[k] = acc_en[1];
				k = k + 1;
			}
			else if(ret_const == 3) //Both children 1 and 2 are accepted.
			{
				if(k == (tot - 1)) //this is the scenario when only one config and its energy can be accpted in population, but there are two accepted children.
				{
					for(i=0;i<num;i++)
                                        {
                                                for(j=0;j<width;j++)
                                                {
                                                        pos_pre[((k * num) + i)][j] = pos1[i][j];
                                                }
                                                for(j=0;j<width_ele;j++)
                                                {
                                                        atom_pre[((k * num) + i)][j] = atom1[i][j];
                                                }
                                        }
					en[k] = acc_en[0];
					k=k+1;
				}
				else
				{
					//printf("here\n");
					for(i=0;i<num;i++) 
                        		{
                                		for(j=0;j<width;j++)
                                		{
                                        		pos_pre[((k * num) + i)][j] = pos1[i][j];
                                		}
                                		for(j=0;j<width_ele;j++)
                                		{		
                                        		atom_pre[((k * num) + i)][j] = atom1[i][j];
                                		}
                        		}
					for(i=0;i<num;i++)
                        		{
                                		for(j=0;j<width;j++)
                                		{
                                        		pos_pre[(((k+1) * num) + i)][j] = pos2[i][j];
                                		}
                                		for(j=0;j<width_ele;j++)
                                		{
                                        		atom_pre[(((k+1) * num) + i)][j] = atom2[i][j];
                                		}
                        		}
					en[k] = acc_en[0];
					en[k+1] = acc_en[1];
					k=k+2;
				}
				//k = k + 2;
			}
			else if(ret_const == 4) //Both children 1 and 2 are rejected.
			{
				k = k + 0;
			}
		} while(k < tot);//TODO Change to tot here. Swapping between two parents starts here. 
		//swapping between two parents here.
                //printf("%lf\t%lf\n",en[0],en[1]);
		//copying configurations, which will be used as parent in next generation.
		//INTRODUCE PARENT SELECTION HERE
		for(i=0;i<tot;i++)
		{
			res[i] = -1;
		}
		//Writing energy values to the ene.$gen file
		//writting energy values 
		char store[15];
		sprintf(store,"ene.%d.csv",gen);

		FILE *fe;
		fe = fopen(store,"a");
		if(fe == NULL)
                {
                        printf("Error in opening file-ene.*\n");
                        exit(1);
                }
		for(i=0;i<tot;i++)
		{
			fprintf(fe,"%lf\n",en[i]);
		}
		fclose(fe);
	
		par_sel(tot,en,res); //populating res[] array, which contains the identity of configurations, which will act as parent in next generation.
		//generating parent files to be used in next generation.
		for(i=0;i<tot;i++)
        	{
                	u = *(res+i);
                	for(j=0;j<num;j++)
                	{
                        	for(k=0;k<width;k++)
                        	{
                                	pos_inp[((i * num) + j)][k] = pos_pre[((num * u) + j)][k];
                        	}
                	}
        	}
        	for(i=0;i<tot;i++)
        	{
                	u = res[i];
                	for(j=0;j<num;j++)
                	{
                        	for(k=0;k<width_ele;k++)
                        	{
                                	atom_inp[((i * num) + j)][k] = atom_pre[((num * u) + j)][k];
                        	}
                	}
        	}
		
		for(i=0;i<tot;i++)
		{
			u = res[i];
			en_par[i] = en[u]; //might not need this.
		}	
		//writting the status of the GA
		status_write(gen);

		//writing checkpoint file
		if(gen%check == 0)
		checkpoint(gen,pos_inp,atom_inp);
	}//ga cycle ends here.
    
	//freeing the memory	
	for(i=0;i<total;i++)
	{
		free(pos[i]);
	}
	for(i=0;i<3;i++)
        {
                free(atom[i]);
        }
	for(i=0;i<total;i++)
        {
                free(pos_inp[i]);
        }
        for(i=0;i<3;i++)
        {
                free(atom_inp[i]);
        }
	
		
exit(0);
}
/*********************************************************************/
//this function generates a random number between interval min_num and max_num, while avoiding a number (ear_v), if it lies in that interval.
int randominlimits(int min_num,int max_num,int ear_v)
{
	int rn = -1;
	int rn1 = -1;

	
	if(ear_v >= min_num && ear_v <= max_num)
	{
		if(ear_v == min_num) //if forbidden number is minimum number 
		{
			rn = rand() % (max_num - min_num - 1) + (min_num + 1);
		}
		else if(ear_v == max_num) //if forbidden number id maximum number.
		{
			rn = rand() % (max_num-min_num-1) + min_num;
		}
		else 
		{
			rn = rand() % (max_num - min_num) + min_num;
			if(rn == ear_v) //if random number is forbidden number
			{	
				rn1 = rand() % 2; //a random number, which is either 0 and 1 is generated.
				if(rn1 == 0)
				{
					rn--; //if random number generated above is 0, then rn is decreased by one.
				}
				else if(rn1 == 1)
				{
					rn++; //if random number generated above is 0, then rn is increased by one.
				}
				else
				{
					;
				}
			}
		
		}
	}
	else //in case forbidden number is not in the specified interval.
	{	
		rn = rand() % (max_num - min_num) + min_num;
	}

	if(rn == ear_v) //just a check!
	{
		printf("Same v as ear_v for %d\n",rn);
		exit(1);
	}
	else
	{
		return rn;
	}
}	
/*********************************************************************/
void parent_gen(int res[])
{
	//TEST
	if(tot < pop_div)
	{
		printf("total number of parent population is less than population division\n");
		exit(1);
	}
	if(tot%pop_div != 0)
	{
		printf("tot is not divisible by pop_div\n");
		exit(1);
	}	
	int a,b,u,v,w,max_num,each_num,counter,x;
	each_num=tot/pop_div;
	//two-dimensional array to store the identity of parents
        int *parent[pop_div];
        for(a=0;a<pop_div;a++)
	{
		parent[a]=(int *)malloc(each_num*sizeof(int));
                if(parent[a] == NULL)
                {       
                        printf("Issue with array allocation- 7\n");
                        exit(1);
                }
	}
	counter=0;
        for(a=0;a<pop_div;a++)
	{
		for(b=0;b<each_num;b++)
		{
			parent[a][b]=counter;
			counter++;
		}
	}
	max_num=(int)pow(div_fact,pop_div);
	//int countnum[tot]; //counting the number of times a particular composition is chosen.
        int ear_v = -1; //variable to store the value of v in particular case, which would be compared with v in next random generation.
        for(a=0;a<tot;a++) 
	{
		u = (rand() % max_num)+1; //finding random number b/n 0 and div_fact^pop_div
		//u = (rand() / (RAND_MAX/max_num + 1)) + 1;
		w=0;
		b=1;
		v=0;
		while(v > 0 || v == 0)
		{
                	v = u - ( w + (max_num/((int)pow(div_fact,b))));
			if(v < 0 || v == 1 || v == 0)
		        {
				x = randominlimits(0,each_num,ear_v);
				res[a] = parent[b-1][x]; //because b starts with 1.
			        ear_v = x;
				break;
			}
			else if( v > 0)
			{
				w = w + (max_num/((int)pow(div_fact,b)));
			}
			b++;
		}
	}
	//freeing the memory
        for(a=0;a<pop_div;a++)
        {
                free(parent[a]);
        }
}
/*****************************************************************/
int swap(char pos1[][width],char pos2[][width],char atom1[][width_ele],char atom2[][width_ele],int gen,double acc_en[],double par_en[])
{
	int i,j,k,l,m,rand_num,num_swap,atom_num,u,n;
	double frac = 0;
	char tempcord1[width];
        char tempat1[width_ele];
        char tempcord2[width];
        char tempat2[width_ele];
	double p = 0.0;
	double q = 0.0;
	int ret_const = -1; // this variable would return a constant, whose value can be 1,2,3 and 4 depicting child 1, child 2, both children and none of children, repsectively have been chosen.
	double delta_ene = 0.0;
	int yon[2];
	atom_num = num;
        rand_num = (rand() % 3 + 1); //the percentage of swapping is kept between 10-30% of total number of atoms in CONFIG file.
        frac = rand_num/10.0;
        num_swap = (int)(atom_num * frac); 
        
	for(i=0;i<num_swap;i++) //change to num_swap here.
        {
                j = rand() % atom_num + 1; //finding 
                for(k=0;k<width;k++)
                {
                        tempcord1[k] = pos1[j][k];
                        tempcord2[k] = pos2[j][k];
                }
                for(k=0;k<width_ele;k++)
                {
                        tempat1[k] = atom1[j][k];
                        tempat2[k] = atom2[j][k];
                }
                //NOTE THAT AT SAME VALUE OF j THE IDENTITY OF ATOM IN BOTH ATOM1[][] AND ATOM2[][] WOULD BE SAME, DUE TO THE WAY CONFIG FILES ARE WRITTEN IN DL_POLY.  
                //Finding the atom number in config.2, which has same coordinates as jth atom in config.1 
                for(l=0;l<atom_num;l++) //change to atom_num here.
                {
                        m = 0;
                        for(n=0;n<width;n++)
                        {
                                if(tempcord1[n] == pos2[l][n])
                                {
                                        m = m+1;
                                }
                        }
                        if(m == width)
                        {
                                break;
                        }
                        else
                        {
                                ;
                        }
                }
                //lth atom in config.2 has same coordinate as jth atom in config.1
		for(u=0;u<atom_num;u++)
                {
                        m = 0;
                        for(n=0;n<width;n++)
                        {
                                if(tempcord2[n] == pos1[u][n])
                                {
                                        m = m+1;
                                }
                        }
                        if(m == width)
                        {
                                break;
                        }
                        else
                        {
                                ;
                        }
                }
                //uth atom in config.1 has same coordinates as jth atom in config.2

                //swapping coordinates of jth and uth atom in config.1 file and jth and uth atom in config.2 file.
                for(k=0;k<width;k++)
                {
                        pos1[j][k] = tempcord2[k];
                        pos1[u][k] = tempcord1[k];

                        pos2[j][k] = tempcord1[k];
                        pos2[l][k] = tempcord2[k];
                }
        }
	
	FILE *f_conf1;
        f_conf1 = fopen("conf.1","w");
        FILE *f_conf2;
        f_conf2 = fopen("conf.2","w");
        
	for(k=0;k<num;k++)
	{
		fprintf(f_conf1,"%d ",k+1); //since the index starts from 1 in the lammps data file
		fprintf(f_conf2,"%d ",k+1);
		for(j=0;j<width_ele-1;j++)
		{ 
			fprintf(f_conf1,"%c",atom1[k][j]);
			fprintf(f_conf2,"%c",atom2[k][j]);
		}
		//introducing space between atom type and coordinates
		fprintf(f_conf1," ");
		fprintf(f_conf2," ");
		for(j=0;j<width;j++)
		{
			fprintf(f_conf1,"%c",pos1[k][j]);
			fprintf(f_conf2,"%c",pos2[k][j]);
		}
	}
        fclose(f_conf1);
        fclose(f_conf2);
	
	//generating two config files
	system("./join.bash");
	FILE *x;
	x = fopen("energy-ga.txt","r");
	if(x == NULL)
	{
		printf("Issue with energy-ga.dat file opening\n");
		exit(1);
	}
	fscanf(x,"%lf",&p);
	fscanf(x,"%lf",&q);
	fclose(x);
	//INSERT ACCEPTANCE CRITERIA HERE
	//yon[] stores the information, whether children are accepted or not. yon[i]=0 or 1 depending whether ith childred is rejectd or accpeted, respectively.
	yon[0] = -1;//initialisation
	yon[1] = -1;
	//checking child acceptance
	delta_ene = p - par_en[0];
	if(delta_ene <= 0.0)
	{
		yon[0] = 1;
	}
	else if(delta_ene > 0.0)
	{
		yon[0] = child_accept(delta_ene);	
	}
	else
	{
		;
	}
	//checking child-2 acceptance
	delta_ene = q - par_en[1];
	if(delta_ene <= 0.0)
        {
                yon[1] = 1;
        }
        else if(delta_ene > 0.0)
        {
		yon[1] = child_accept(delta_ene);
        }
        else
        {
                ;
        }
        //determining value of ret_const
	if(yon[0] == 1 && yon[1] == 0)
	{
		ret_const = 1;
	}
	else if(yon[0] == 0 && yon[1] == 1)
	{
		ret_const = 2;
	}
	else if(yon[0] == 1 && yon[1] == 1)
	{
		ret_const = 3;
	}
	else if(yon[0] == 0 && yon[1] == 0)
        {
                ret_const = 4;
        }
	else
	{
		printf("Issue with ret_const determination\n");
		exit(1);
	}
	/*****************************************/	
	//ret_const = 3;  //The value of 3 has been simply assigned, in case of acceptance criteria has not been defined and all the children are accepted.

	//this section should be invoked if only accepted energy values need to be added.
	if(ret_const == 1)
	{
		acc_en[0] = p;	
	}
	else if(ret_const == 2)
	{
		acc_en[1] = q;
	}
	else if(ret_const == 3)
	{
		acc_en[0] = p;
		acc_en[1] = q;
	}
	else if(ret_const == 4)
	{
		;
	}
	else
	{
		;
	}
	
return ret_const;
}
/*************************************************************************************************/
void status_write(int fil_num)
{
	int i;
        int lin_num = -1;
        int degen_num = 0;
        char name[15];
        double mean_ene = 0.0;
        double sd_ene = 0.0;
        
	FILE *f_ga;
        f_ga = fopen("status_ga.csv","a");
        if(f_ga == NULL)
        {
                printf("Issue with status_ga.csv file opening\n");
                exit(1);
        }
	sprintf(name,"ene.%d.csv",fil_num);
        //count number of lines in file
        lin_num = countlines(name);
        

        double *check_ene;
        check_ene = malloc(lin_num * sizeof(double));
        if(check_ene == NULL)
        {
                printf("Issue with array allocation- 1\n");
                exit(1);
        }
        FILE *fene;
        fene = fopen(name,"r");
        if(fene == NULL)
        {
                printf("Issue with ene.* file opening\n");
                exit(1);
        }
        for(i=0;i<lin_num;i++)
        {
                fscanf(fene,"%lf",&check_ene[i]);
        }
        fclose(fene);
        //mean energy claculation
        mean_ene = 0.0;
        for(i=0;i<lin_num;i++)
        {
                mean_ene = mean_ene + check_ene[i];
        }
        mean_ene = mean_ene/lin_num;
        //standard deviation calculation
        sd_ene = 0.0;
        for(i=0;i<lin_num;i++)
        {
                sd_ene = sd_ene + ((mean_ene - check_ene[i])*(mean_ene - check_ene[i]));
        }
        sd_ene = sd_ene/(lin_num - 1);
        sd_ene = sqrt(sd_ene);

        ascend(check_ene,lin_num);
        degen_num = degen_count(check_ene,lin_num);
 	
	fprintf(f_ga,"%d\t%lf\t%lf\n",degen_num,mean_ene,sd_ene);       
	
	fclose(f_ga);
	free(check_ene);
}
/*************************************************************************************************/
int countlines(char *filename)
{
        int ch = 0;
        int lines = 0;
        FILE *fp;
        fp = fopen(filename,"r");
        if(fp == NULL)
        {
                printf("Issue with opening file for counting the number of lines\n");
                exit(0);
        }
        //lines++;
        while ((ch = fgetc(fp)) != EOF)
        {
                if (ch == '\n')
                        lines++;
        }
        fclose(fp);
  	return lines;

}
/************************************************************************************************/
void ascend(double check_ene[],int lin)
{
        int i,j;
        double a = 0.0;

        for(i=0;i<lin;i++)
        {
                for(j=i+1;j<lin;j++)
                {
                        if(check_ene[i] > check_ene[j])
                        {
                                a = check_ene[i];
                                check_ene[i] = check_ene[j];
                                check_ene[j] = a;
                        }
                }
        }
}
/*************************************************************************************************/
int degen_count(double check_ene[],int lin)
{
        int i,count;

        count = 0;
        for(i=0;i<lin-1;i++)
        {
                if(check_ene[i] != check_ene[i+1])
		{
                	count++;
		}
        }
	return count;
}
/**************************************************************************************************/
void par_sel(int totl,double check_ene[],int res[])
{
	int i,j,k,l,yon;
	double b = 0.0;
	double Q = 0.0;
	double beta = 1/(rb*T);
	double *ene;
	ene = malloc(totl*sizeof(double));
	if(ene == NULL)
	{
		printf("Issue with array allocation in par_sel function\n");
		exit(1);
	}
	double *prob;
	prob = malloc(totl*sizeof(double));
        if(prob == NULL)
        {
                printf("Issue with array allocation in par_sel function\n");
                exit(1);
        }	

	for(i=0;i<totl;i++)
        {
                ene[i] = check_ene[i]/num;
	}
	//finding minimum energy
        for(i=0;i<totl;i++)
        {
                for(j=i+1;j<totl;j++)
                {
                        if(ene[i] > ene[j])
                        {
                                b= ene[i];
                                ene[i] = ene[j];
                                ene[j] = b;
                        }
                }
        }
        for(i=0;i<totl;i++)
        {
                ene[i] = ene[i]*96.845;
        }

        for(i=0;i<totl;i++)
        {
                ene[i] = ene[i]*(-1)*beta;
        }

        for(i=0;i<totl;i++)
        {
                ene[i] = exp(ene[i]);
	}
        Q = 0.0;
        for(i=0;i<totl;i++)
        {
                Q = Q + ene[i];
        }
        for(i=0;i<totl;i++)
        {
                prob[i] = ene[i]/Q;
        }
	
	//filling up res[] array
        i = 0;
        while(i < tot)//change to tot here.
        {
                rep: k = rand()%totl;
                yon = chooseornot(k,prob);
                if(yon == 0)
                {
                        goto rep;
                }

                yon = -1;
                rep1: l = randominlimits(0,totl,k);
                yon = chooseornot(l,prob);
                if(yon == 0)
                {
                        goto rep1;
                }
                res[i] = k;
                res[i+1] = l;
                i = i +2;
        }
	
free(ene);
free(prob);
}
/*************************************************************************************************/
int chooseornot(int u,double prob[])
{
        double t = 0.0;
        int fact = 1000000;
        int tar = 0;
        int yon = -1;
        t = prob[u];
        int rn = -1;
        tar = fact*t;
        rn = rand()%fact+1;
        if(rn >= 1 && rn <= tar)
        {
                yon = 1;
        }
        else
        {
                yon = 0;
        }
        return yon;
}
/*************************************************************************************************/
int child_accept(double delta_ene)
{
	double beta = 1/(rb*T);
	double ene = 0.0;
	double pr = 0.0;
	int fact = 1000000;
	int tar = 0;
	int yn = -1;
	int rn = -1;

	ene = delta_ene/num;
	ene = ene*96.845;
	
	pr = exp((-1)*ene*beta);
	tar = fact*pr;
        rn = rand()%fact+1;
        if(rn >= 1 && rn <= tar)
        {
                yn = 1;       
        }
        else
        {
                yn = 0;
        }
        return yn;
	
}
/*************************************************************************************************/
void checkpoint(int gen,char *pos_inp[],char *atom_inp[])
{
	int a,b,c;
	system("rm -f alist.backup.*");
	system("rm -f coord.backup.*");
	char name1[20]; //keeping the name of the file within 15 char
	char name2[20];
	for(a=0;a<tot;a++)
	{
		sprintf(name1,"cord.backup.%d",a);
		FILE *fx;
		fx = fopen(name1,"w");
                if(fx == NULL)
                {
                        printf("Error in opening file-backup.*\n");
                        exit(1);
                }
                for(b=0;b<num;b++)
                {
			for(c=0;c<width;c++)
			{
				fprintf(fx,"%c",pos_inp[((a * num) + b)][c]);
			}
                }
                fclose(fx);
		sprintf(name2,"alist.backup.%d",a);
		FILE *fy;
		fy = fopen(name2,"w");
		if(fy == NULL)
		{
			printf("Error in opening file-backup.*\n");
                        exit(1);
		}
		for(b=0;b<num;b++)
                {
                        for(c=0;c<width_ele;c++)
                        {
                                fprintf(fy,"%c",atom_inp[((a * num) + b)][c]);
                        }
                }
		fclose(fy);
	}
}
