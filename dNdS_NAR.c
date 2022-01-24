//-----------------------------Find Syn/NonSyn Mutation, (dN/dS values) by pairwise comparison of two input sequences -------------------------
//Program is written in C programing language and compiled using code blocks version 12.11 in windows 10 operating system
/*Sample gene sequences are given in two input files in fasta format, File1_Ec.txt and File2_Se.txt
File1_Ec.txt contains 32 gene sequences of E.coli and
File2_Se.txt contains 32 gene sequences of S. Enterica
In both the files genes are arranged in same order.
Corresponding genes of both the organisms are equal in size.
Program calculates result as No. of Syn Mutations, NonSyn Mutations, SynSite, NonSynSite and dN/dS
and stores the result in a output file dNdS_Result.txt
*/
#include <stdio.h>
int findCodon(char, char, char);
int main()
{
	FILE *in1, *in2, *out2;
	int i,j,len,mx;
	char c,fNuc,str[100],info1[100][20],info2[100][20],seq[10000], refSeq1[100][5000],refSeq2[100][5000],mutNuc[40000];
	int size[100];

    char AA3[64][4]={"Phe","Ser","Tyr","Cys","Phe","Ser","Tyr","Cys","Leu","Ser","TER","TER","Leu","Ser","TER","Trp","Leu","Pro","His","Arg","Leu","Pro","His","Arg","Leu","Pro","Gln","Arg","Leu","Pro","Gln","Arg","Ile","Thr","Asn","Ser","Ile","Thr","Asn","Ser","Ile","Thr","Lys","Arg","Met","Thr","Lys","Arg","Val","Ala","Asp","Gly","Val","Ala","Asp","Gly","Val","Ala","Glu","Gly","Val","Ala","Glu","Gly"};
    char AA1[64]={'F','S','Y','C','F','S','Y','C','L','S','X','X','L','S','X','W','L','P','H','R','L','P','H','R','L','P','Q','R','L','P','Q','R','I','T','N','S','I','T','N','S','I','T','K','R','M','T','K','R','V','A','D','G','V','A','D','G','V','A','E','G','V','A','E','G'};
    char Codon[64][4]={"UUU","UCU","UAU","UGU","UUC","UCC","UAC","UGC","UUA","UCA","UAA","UGA","UUG","UCG","UAG","UGG","CUU","CCU","CAU","CGU","CUC","CCC","CAC","CGC","CUA","CCA","CAA","CGA","CUG","CCG","CAG","CGG","AUU","ACU","AAU","AGU","AUC","ACC","AAC","AGC","AUA","ACA","AAA","AGA","AUG","ACG","AAG","AGG","GUU","GCU","GAU","GGU","GUC","GCC","GAC","GGC","GUA","GCA","GAA","GGA","GUG","GCG","GAG","GGG"};
    float SynSite[64]={0.666666667,1,0.666666667,0.666666667,0.666666667,1,0.666666667,0.666666667,1.333333333,1,-1,-1,1.333333333,1,-1,0,1,1,0.666666667,1,1,1,0.666666667,1,1.666666667,1,0.666666667,1.166666667,1.666666667,1,0.666666667,1.166666667,0.833333333,1,0.666666667,0.666666667,0.833333333,1,0.666666667,0.666666667,0.333333333,1,0.666666667,0.833333333,0,1,0.666666667,0.833333333,1,1,0.666666667,1,1,1,0.666666667,1,1,1,0.666666667,1,1,1,0.666666667,1};
    float NonSynSite[64]={2.333333333,2,2,2.166666667,2.333333333,2,2,2.166666667,1.333333333,1.666666667,-1,-1,1.5,1.833333333,-1,1.666666667,2,2,2.333333333,2,2,2,2.333333333,2,1.333333333,2,1.666666667,1.166666667,1.333333333,2,1.666666667,1.833333333,2.166666667,2,2.333333333,2.333333333,2.166666667,2,2.333333333,2.333333333,2.666666667,2,2.167,2,3,2,2.167,2.166666667,2,2,2.333333333,2,2,2,2.333333333,2,2,2,2.167,1.833333333,2,2,2.167,2};

    int codeR,codeM,synMut,nonSynMut;
    float synSite1,nonSynSite1,synSite2,nonSynSite2,dNdS;


	int n,k,l;

	//Preprocessing-Read input gene sequences, calculate size
	in1=fopen("File1_Ec.txt","r");
    n=0;
    while((c=getc(in1))!=EOF)
        if(c=='>')n++;
    fclose(in1);
    printf("n=%d\n",n);

    in1=fopen("File1_Ec.txt","r");
	in2=fopen("File2_Se.txt","r");

	out2=fopen("dNdS_Result.txt","w");

	for(k=0;k<n;k++){
		i=0;
		while((c=getc(in1))!='\n'){
            info1[k][i]=c;i++;
		}
		info1[k][i]='\0';
		i=0;while((c=getc(in1))!='\n'){
		    refSeq1[k][i]=c;i++;
		}
		refSeq1[k][i]='\0';
		size[k]=i;
	}
	fclose(in1);

	for(k=0;k<n;k++){
		i=0;
		while((c=getc(in2))!='\n'){
            info2[k][i]=c;i++;
		}
		info2[k][i]='\0';
		i=0;while((c=getc(in2))!='\n'){
		    refSeq2[k][i]=c;i++;
		}
		refSeq2[k][i]='\0';
	}
	fclose(in2);
    //----End Preprocessing--------------
    //Find syn/nonSyn Mutations and calculate dN/dS
    for(k=0;k<n;k++){
        synMut=nonSynMut=0;
        printf("Size:%d\n",size[k]);
        for(i=0;i<size[k];i++){
            if(refSeq1[k][i]!=refSeq2[k][i]){
                if(i%3==0){
                    codeR=findCodon(refSeq1[k][i],refSeq1[k][i+1],refSeq1[k][i+2]);
                    codeM=findCodon(refSeq2[k][i],refSeq1[k][i+1],refSeq1[k][i+2]);
                }
                else if (i%3==1){
                    codeR=findCodon(refSeq1[k][i-1],refSeq1[k][i],refSeq1[k][i+1]);
                    codeM=findCodon(refSeq1[k][i-1],refSeq2[k][i],refSeq1[k][i+1]);
                }
                else if (i%3==2){
                    codeR=findCodon(refSeq1[k][i-2],refSeq1[k][i-1],refSeq1[k][i]);
                    codeM=findCodon(refSeq1[k][i-2],refSeq1[k][i-1],refSeq2[k][i]);
                }
                if(codeR<65&&codeM<65){
                    if(AA1[codeR-1]==AA1[codeM-1])synMut++;
                    else nonSynMut++;
                }
            }
        }
        synSite1=nonSynSite1=synSite2=nonSynSite2=0;
        for(i=0;i<size[k];i=i+3){
            codeR=findCodon(refSeq1[k][i],refSeq1[k][i+1],refSeq1[k][i+2]);
            codeM=findCodon(refSeq2[k][i],refSeq2[k][i+1],refSeq2[k][i+2]);
            if(codeR<65){
                synSite1+=SynSite[codeR-1];
                nonSynSite1+=NonSynSite[codeR-1];
            }
            if(codeM<65){
                synSite2+=SynSite[codeM-1];
                nonSynSite2+=NonSynSite[codeM-1];
            }
        }
        dNdS=(nonSynMut/((nonSynSite1+nonSynSite2)/2.0)/(synMut/((synSite1+synSite2)/2.0)));
        fprintf(out2,"%s\t%d\t%d\t%f\t%f\t%f\n",info1[k],synMut,nonSynMut,(synSite1+synSite2)/2,(nonSynSite1+nonSynSite2)/2,dNdS);
  	}
	fprintf(out2,"\nGene\tNo. Syn Mutation\tNo. NonSyn Mutation\tSynSite\tNonSynSite\tdN/dS\n");
	fclose(out2);
	return(0);
}

int findCodon(char c1, char c2, char c3){
    int pc[3],code;
    pc[0]=pc[0]=pc[0]=100;

    if(c1=='u'||c1=='t'||c1=='U'||c1=='T') {pc[0]=1;}
    else if(c1=='c'||c1=='C'){pc[0]=2;}
    else if(c1=='a'||c1=='A'){pc[0]=3;}
    else if(c1=='g'||c1=='G'){pc[0]=4;}

    if(c2=='u'||c2=='t'||c2=='U'||c2=='T') {pc[1]=1;}
    else if(c2=='c'||c2=='C'){pc[1]=2;}
    else if(c2=='a'||c2=='A'){pc[1]=3;}
    else if(c2=='g'||c2=='G'){pc[1]=4;}

    if(c3=='u'||c3=='t'||c3=='U'||c3=='T') {pc[2]=1;}
    else if(c3=='c'||c3=='C'){pc[2]=2;}
    else if(c3=='a'||c3=='A'){pc[2]=3;}
    else if(c3=='g'||c3=='G'){pc[2]=4;}

    code=(pc[0]-1)*16+pc[1]+(pc[2]-1)*4;
    return code;
}

