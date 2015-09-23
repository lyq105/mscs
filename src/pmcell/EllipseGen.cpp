#include "EllipseGen.h"
#define pi 3.1415926
#define N1 (int)(2*(upbar-downbar)/(b_aver))
#define N2 (int)(2*(rightbar-leftbar)/(b_aver))
#define N3 N1*N2
#define h 0.05


using namespace std;

	struct parameter0
{
	double al,bl,cl,a,b,c,center_x,center_y,center_z,ymin,ymax,xmin,xmax,zmin,zmax;
	double alpha1,bata1,gama1,alpha2,bata2,gama2,alpha3,bata3,gama3;
	double gamamin,gamamax,miumin,miumax;
	double A,B,C;
	double coe[7];
	
};

struct parameter
{
	double al,bl,cl,a,b,c,center_x,center_y,center_z;
	double A,B,C;
	double alpha1,bata1,gama1,alpha2,bata2,gama2,alpha3,bata3,gama3;
	struct parameter *next;
};

struct basicparameter
{
	double al,bl,cl,a,b,c,center_x,center_y,center_z;
	double alpha1,bata1,gama1,alpha2,bata2,gama2,alpha3,bata3,gama3;
	double A,B,C;
};

struct coordinate
{
	double x;
	double y;
	double z;
};

char percent='%';                      
double M=pow(2.0,42);
double lamda=pow(5.0,17);
struct parameter *head;
double b_aver=0,a_aver=0;
double famin=0,famax=0,fb=0;
double eamin=0,eamax=0,ebmin=0,ecmin=0;
double camin=0,camax=0,cbmin=0,cc=0;
double gamamin=0,gamamax=0,miumin=0,miumax=0;
double sitagama=0,miugama=0,sitamiu=0,miumiu=0;
double leftbar,rightbar,formbar,backbar,upbar,downbar,leftbar1,rightbar1,backbar1,formbar1,downbar1,upbar1;
int initialx,finalx,initialy[2],finaly[2];
int initialwz[2],finalwz[2],initialkz[2],finalkz[2];
double *total_xmin[2],*total_xmax[2],*total_ymin[2],*total_ymax[2];
double error=pow(10.0,-8),errort=pow(10.0,-8);
double f_xseed,f_yseed,f_zseed,f_aseed;
double e_xseed,e_yseed,e_zseed,e_aseed,e_bseed,e_cseed;
double c_xseed,c_yseed,c_zseed,c_aseed,c_bseed;
double 	ig1seed,im1seed,ig2seed,im2seed;
double mseed;
double seqx,seqy,seqz;


int EllipseGen::uniell_generation(string input_file, string out_file, string data_file)
{	
	srand((unsigned int)time(0));
	f_xseed=seed();  f_yseed=seed();  f_zseed=seed();  f_aseed=seed();
	e_xseed=seed();  e_yseed=seed();  e_zseed=seed();  e_aseed=seed();  e_bseed=seed();  e_cseed=seed();
	c_xseed=seed();  c_yseed=seed();  c_zseed=seed();  c_aseed=seed();  c_bseed=seed();
	ig1seed=seed();  im1seed=seed();  ig2seed=seed();  im2seed=seed();
	mseed=seed();
	seqx=seed();  seqy=seed();  seqz=seed();


	char rline[250];
	int distribution_sign1=0,distribution_sign2=0;
	ifstream input_stream;
	
	input_stream.open( input_file.c_str(), ios::in );

	if(!input_stream)
	{
        cout << "Can not open the file: "  << "\n";
        return 0 ;          
     };

	istrstream istr0(get_line(input_stream,rline));
        istr0 >> leftbar >> rightbar >> backbar >> formbar
			  >> downbar >> upbar;//长方体区域的六条边
	

	istrstream istr2(get_line(input_stream,rline));
        istr2 >> distribution_sign1;//倾角分布信息，取0为均匀分布，取1为正态分布


	//ofstream output_stream;
	//output_stream.open(out_file.c_str(),ios::out);

 //   if(!output_stream)
	//{
	//	cout <<"Cannot open file.\n" <<out_file <<endl;
	//	return 0;
	//}

	//output_stream <<"%椭球的个数" <<endl;
	//output_stream <<1<<endl;

 //   //输出单胞位置和尺寸
	//output_stream <<"%单胞六个边" <<endl;
	//output_stream  <<leftbar <<" " <<rightbar<<" "<<backbar<<" "<<formbar <<" "<<downbar <<" "
	//	          <<upbar <<endl;

 ////输出椭球18个参数（参数方程：椭球中心坐标（3个）、边界层厚度（3个）、尺寸（3个）、方位角（9个））
	//output_stream <<"%输出椭球18个参数:椭球中心坐标、边界层厚度、尺寸、方位角" <<endl;

	//for(int k=0;k<1;k++)
	//{
	//	for(int j=0;j<1;j++)
	//	{
	//		for(int i=0;i<1;i++)
	//		{
	//			double c_x,c_y,c_z;
	//			c_x=0.5+i*0.5;
	//			c_y=0.5+j*0.5;
	//			c_z=0.5+k*0.5;
 //                output_stream  <<c_x <<" " <<c_y <<" " <<c_z <<" "
	//			<<0.47 <<" " <<0.47 <<" " <<0.47<<" "<<0.47 <<" " <<0.47 <<" " <<0.47<<" "
	//			<<1 <<" " <<0 <<" " <<0 <<" "
 //               <<0 <<" " <<1 <<" " <<0 <<" "
	//			<<0 <<" " <<0 <<" " <<1 <<endl;
	//		}
	//	}
	//}


	//neededrate=4*pi*pow(0.4705,3)/3;
	//realrate=neededrate;
	//output_stream.close();

	compact_elli(input_stream,distribution_sign1);

	
	istrstream istr3(get_line(input_stream,rline));
        istr3 >> distribution_sign2;//中心分布信息，取0为均匀分布，1为正态分布，2为指数分布

	/*if(distribution_sign2==0)
    	uniform(input_stream,out_file);
	else*/
	if(no_uniform(distribution_sign2,input_stream,out_file, data_file)==0)
	{
		hout << "可能是体积分数不够！" << endl;
		return 0;
	}
	
	input_stream.close();

	return 1;
};

int EllipseGen::compact_elli(ifstream& input_stream,int distribution_sign1)
{
	
	int i,j,k,k1;
	int i1,i2,j1;
	int fm,em,cm,m;
	int e_num=0,f_num=0,c_num=0;
	int n,num=0;
	int norm1,norm2,norm4,norm5;
	int f_e_c_division;    //若要生成纤维取0，若要生成椭球取1，若要生成薄钱币颗粒取2
	double *ixseed,*iyseed,*izseed;
	double dz,dy,dx;
	char rline[250];
	int shape_sign1,shape_sign2,shape_sign3;//shape1取1表示存在纤维颗粒，取0表示不存在
                                                // 同理，shape2,shape3取1代表椭球颗粒以及薄钱币颗粒存在，否则不存在
	struct parameter0 **elliptic[3],*p0;
	struct parameter *p1,*p2;

/***********************确定有几种形态的颗粒，并读取相关的参数*************************/

   istrstream istr1(get_line(input_stream,rline));
   istr1 >>shape_sign1 >>shape_sign2 >>shape_sign3 ;
 
  /*输入纤维的amin amax bmin cmin gamamin gamamax miumin miumax*/
//纤维长轴a服从(famin,famax)的均匀分布,中轴b取fb,短轴c与b同
 //长轴a与z轴的夹角服从(fgamamin,fgamamax)的均匀分布,a在xoy面上投影与x轴的夹角服从(fmiumin,fmiumax)的均匀分布 
  if(shape_sign1==1)
  {
	  istrstream istr20(get_line(input_stream,rline));
      istr20 >> famin >>famax >>fb ;
  }
  /*输入椭球体的amin amax bmin cmin gamamin gamamax miumin miumax*/
//椭球长轴a服从(eamin,eamax)的均匀分布,中轴b服从(ebmin,a)的均匀分布,短轴c服从(ecmin,b)的均匀分布
 //长轴a与z轴的夹角服从(egamamin,egamamax)的均匀分布,a在xoy面上投影与x轴的夹角服从(emiumin,emiumax)的均匀分布
  if(shape_sign2==1)
  {
	  istrstream istr21(get_line(input_stream,rline));
      istr21 >> eamin >>eamax >>ebmin >>ecmin ;
	  
  }
  /****************************************/
  /*输入薄钱币颗粒的amin amax bmin cmin gamamin gamamax miumin miumax*/
//薄钱币颗粒长轴a服从(camin,camax)的均匀分布,中轴b服从(cbmin,a),短轴c取cc
 //长轴a与z轴的夹角服从(cgamamin,cgamamax)的均匀分布,a在xoy面上投影与x轴的夹角服从(cmiumin,cmiumax)的均匀分布
  if(shape_sign3==1)
  {
	  istrstream istr22(get_line(input_stream,rline));
      istr22 >> camin >>camax >>cbmin >>cc;
  }

/*********************颗粒形态确定完毕**********************************************/

 //颗粒方向参数的读取
  if(distribution_sign1==0)
  {
	  istrstream istr30(get_line(input_stream,rline));
      istr30 >>gamamin >>gamamax >>miumin >>miumax;
  }
  else
  {
	  istrstream istr31(get_line(input_stream,rline));
      istr31 >>sitagama >>miugama >>sitamiu >>miumiu;

  }

  //纤维，椭球，薄钱币颗粒的个数比为fm:em:cm
  istrstream istr4(get_line(input_stream,rline));
      istr4 >> fm >>em >>cm ;
  /*****************输入长方体区域参数**********************************/
  
 
  m=fm+em+cm;
  
  
  b_aver=(fb*fm/(fm+em+cm)+ebmin*em/(fm+em+cm)+cbmin*cm/(fm+em+cm)); //每一行可能容纳颗粒的最多个数由其决定 

  int sum=0;
  if(shape_sign1==1)
  {
	  sum=sum+3;
	  a_aver=a_aver+(famin+famax+fb);
  }
  if(shape_sign2==1)
  {
	  sum=sum+3;
	  a_aver=a_aver+(eamin+eamax+ebmin);
  }
  if(shape_sign3==1)
  {
	  sum=sum+3;
	  a_aver=a_aver+(camin+camax+cbmin);
  }
  
  a_aver=a_aver/sum;
  
 
  

  leftbar1=leftbar;
  rightbar1=rightbar;
  downbar1=downbar;
  upbar1=upbar;
  backbar1=backbar;
  formbar1=formbar;
  //申请空间

//elliptic[2]存储待生成层的颗粒，elliptic[0]，elliptic[1]存储此层前1层，前2层的颗粒

  for(k=0;k<=2;k++)
  {
      elliptic[k]=(struct parameter0 **)malloc(N1*sizeof(struct parameter0 *));
	  if(elliptic[k]==NULL)
	  {
    	  printf("Can't obtain the needed space");
    	  return 0;
	  }
     for(j=0;j<N1;j++)
	 {
       	  elliptic[k][j]=(struct parameter0 *)malloc(N2*sizeof(struct parameter0));
    	  if(elliptic[k][j]==NULL)
		  {
        	  printf("Can't obtain the needed space");
        	  return 0;
		  }
	 }
  }

  for(k1=0;k1<=1;k1++)
  {
     total_xmin[k1]=(double *)malloc(N2*sizeof(double));//存储待生成层的前1层，前2层每列x的最小值
     total_xmax[k1]=(double *)malloc(N2*sizeof(double));//存储待生成层的前1层，前2层每列x的最大值
     total_ymin[k1]=(double *)malloc(N1*sizeof(double));//存储待生成层的前1层，前2层每列y的最小值
     total_ymax[k1]=(double *)malloc(N1*sizeof(double));//存储待生成层的前1层，前2层每列y的最大值
  }
  //清零
  for(k1=0;k1<=1;k1++)
	  for(i1=0;i1<N2;i1++)
	  {
		  total_xmin[k1][i1]=0;
		  total_xmax[k1][i1]=0;
	  }
  for(k1=0;k1<=1;k1++)
	  for(j1=0;j1<N1;j1++)
	  {
		  total_ymin[k1][j1]=0;
		  total_ymax[k1][j1]=0;
	  }


  for(k=0;k<=2;k++)
	  for(j=0;j<N1;j++)
		  for(i=0;i<N2;i++)
		  {
			  elliptic[k][j][i].A=0;
			  elliptic[k][j][i].a=0;
			  elliptic[k][j][i].al=0;
			  elliptic[k][j][i].alpha1=0;
			  elliptic[k][j][i].alpha2=0;
			  elliptic[k][j][i].alpha3=0;
			  elliptic[k][j][i].B=0;
			  elliptic[k][j][i].b=0;
			  elliptic[k][j][i].bata1=0;
			  elliptic[k][j][i].bata2=0;
			  elliptic[k][j][i].bata3=0;
			  elliptic[k][j][i].bl=0;
			  elliptic[k][j][i].C=0;
			  elliptic[k][j][i].c=0;
			  elliptic[k][j][i].center_x=0;
			  elliptic[k][j][i].center_y=0;
			  elliptic[k][j][i].center_z=0;
			  elliptic[k][j][i].cl=0;

			  for(k1=0;k1<=6;k1++)
			  {
	    		  elliptic[k][j][i].coe[k1]=0;
			  }

			  elliptic[k][j][i].gama1=0;
			  elliptic[k][j][i].gama2=0;
			  elliptic[k][j][i].gama3=0;
			  elliptic[k][j][i].xmax=0;
			  elliptic[k][j][i].ymax=0;
			  elliptic[k][j][i].xmin=0;
			  elliptic[k][j][i].ymin=0;
			  elliptic[k][j][i].zmin=0;
			  elliptic[k][j][i].zmax=0;
		  }
  

	 n=0; 
	 head=NULL;
//////****************************************** 椭球的定位******************************************************//////
     for(k=0;;k++)
	 {
		 
		 norm1=1;//颗粒的z最大值大于upbar1时 为1，用于判断椭球生成是否完毕
		
		for(j=0;;j++)
		{
			 norm2=1;//颗粒y最大值大于formbar1时为1，用于判断此层椭球是否生成完毕
			for(i=0;;i++)
			{
				norm4=1;//若椭球满足插入下一层的条件则取0
				norm5=1;//若椭球满足插入下一行的条件则取0
				//形成链表
				{
					n=n+1;
					if(n==1)
					{
						p1=p2=(struct parameter *)malloc(sizeof(struct parameter));
				    	head=p1;
					}
					else
					{
						p1=(struct parameter *)malloc(sizeof(struct parameter));
						(*p2).next=p1;
			    		p2=p1;	
					}
				}



				p0=&elliptic[2][j][i];

			

				double m1;
			
				m1=unifrnd(0,(double)(m),&mseed);
	
			
		    	if(m1<fm)
				{
		    	  	f_e_c_division=0;
		    		f_num++;
				}
	    		else if(m1>=fm&&m1<(em+fm))
				{
		     		f_e_c_division=1;
			   		e_num++;
				}
				else
				{
					f_e_c_division=2;
		    		c_num++;
				}
				//printf("1j=%d i=%d N1=%d N2=%d\n",j,i,N1,N2);
           	  parametergeneration1(p0,f_e_c_division,distribution_sign1);	
			  

				parametergeneration2(p0);

				if(f_e_c_division==0)
				{
					ixseed=&f_xseed;
					iyseed=&f_yseed;
					izseed=&f_zseed;
				}
				else if(f_e_c_division==1)
				{
					ixseed=&e_xseed;
					iyseed=&e_yseed;
					izseed=&e_zseed;
				}
				else
				{
					ixseed=&c_xseed;
					iyseed=&c_yseed;
					izseed=&c_zseed;
				}


///*********************************************x的第一次定位******************************************************///
				if(i==0)
				{//printf("%lf %lf %lf %lf\n",cos((*p0).alpha1),cos((*p0).alpha2),(*p0).alpha3,(*p0).A);
				
					(*p0).center_x=unifrnd(leftbar1,leftbar1+a_aver,ixseed);
                    (*p0).xmin=(*p0).center_x-sqrt((*p0).A);
					(*p0).xmax=(*p0).center_x+sqrt((*p0).A);
					
				}
				else
				{
					for(i1=0;i1<i;i1++)
						if(elliptic[2][j][i1].xmax>elliptic[2][j][i-1].xmin)
						{
							initialx=i1;
							break;
						}
					finalx=i-1;//与当前椭球定位相关的当前行椭球，从第initialx个到第finalx个

					(*p0).center_x=elliptic[2][j][i-1].xmax+sqrt((*p0).A);
					(*p0).xmin=(*p0).center_x-sqrt((*p0).A);
					(*p0).xmax=(*p0).center_x+sqrt((*p0).A);
						         
				}
///*********************************************x的第一次定位结束******************************************************///

////*********************************************y的第一次定位******************************************************///
				if(j==0)
				{
					(*p0).center_y=unifrnd(backbar1,backbar1+a_aver,iyseed);
                    (*p0).ymin=(*p0).center_y-sqrt((*p0).B);
					(*p0).ymax=(*p0).center_y+sqrt((*p0).B);
				}
				else
				{
                    //j为当前行，与当前颗粒定位相关的第j-2行颗粒，从第initialy[0]个到第finaly[0]个
					if(j>=2)
					{
						if(i==0)
							initialy[0]=0;
						else
						{
				    		for(i1=0;;i1++)
			    		     	if(elliptic[2][j-2][i1].xmax>elliptic[2][j][i-1].xmin)
								{
			    	    	 		initialy[0]=i1;
			    	    	 		break;
								}
						}

                        for(i1=initialy[0];;i1++)
						{
					        if(elliptic[2][j-2][i1].xmin>elliptic[2][j][i].xmax)
							{
								if(i1==0)
									finaly[0]=0;
								else
						        	finaly[0]=i1-1;
					    		break;
							}
			    			else if(elliptic[2][j-2][i1+1].al==0)
							{
		    					finaly[0]=i1;
		    					break;
							}
						}
					}

                   if(j>=1)
					{////j为当前行，与当前颗粒定位相关的第j-1行颗粒，从第initialy[1]个到第finaly[1]个
		    			if(i==0)
		    				initialy[1]=0;
		    			else
						{
		    	    		for(i2=0;;i2++)
		    	     		   	if(elliptic[2][j-1][i2].xmax>elliptic[2][j][i-1].xmin)
								{
			        	    		initialy[1]=i2;
			        	    		break;
								}
						}
                        for(i2=initialy[1];;i2++)
						{
		    			    if(elliptic[2][j-1][i2].xmin>elliptic[2][j][i].xmax)
							{
		    					if(i2==0)
			    					finaly[1]=0;
			    				else
		 	    		        	finaly[1]=i2-1;
		    		    		break;
							}
			     		    else if(elliptic[2][j-1][i2+1].al==0)
							{
		        				finaly[1]=i2;
		        				break;
							}
						}
					}
					firstsettley(elliptic[2],j,i);
                    (*p0).ymin=(*p0).center_y-sqrt((*p0).B);
					(*p0).ymax=(*p0).center_y+sqrt((*p0).B);
				}
///********************************************y的第一次定位结束*********************************************///
				
///***********************************************z的第一次定位************************************************///

/*********************以下确定当前椭球所能移动的范围，x0,x1为x范围，y0,y1为y的范围*******************************************/
				if(k==0)
				{
					(*p0).center_z=unifrnd(downbar1,downbar1+a_aver,izseed);
					//printf("%d %d %d %lf %lf %lf\n",k,j,i,(*p0).center_z,(*p0).center_y,(*p0).center_x);
					(*p0).zmin=(*p0).center_z-sqrt((*p0).C);
					(*p0).zmax=(*p0).center_z+sqrt((*p0).C);
				}
				else
				{
					double x0,y0,x1,y1;
				
					if(i>=1)
			    		x0=elliptic[2][j][i-1].xmin;
					else
                        x0=elliptic[2][j][i].xmin;

					x1=elliptic[2][j][i].xmax;
					y1=elliptic[2][j][i].ymax;

					
					if(i>=1)
				        y0=min(elliptic[2][j][i-1].ymin,elliptic[2][j][i].ymin);
					else
						y0=elliptic[2][j][i].ymin;

					if(j>=1)
					{
                        for(i1=initialy[1];i1<=finaly[1];i1++)
							if(elliptic[2][j-1][i1].ymin<y0)
								y0=elliptic[2][j-1][i1].ymin;
					}
					/*if(j>=2)
					{
                        for(i1=initialy[0];i1<=finaly[0];i1++)
							if(elliptic[2][j-2][i1].ymin<y0)
								y0=elliptic[2][j-2][i1].ymin;
					}*/

/*********************以下确定第k-1行椭圆中与此椭圆定位有关的所有椭圆*******************************************/
			    	for(i1=0;i1<N2;i1++)
					{
		   				if(total_xmax[1][i1]>x0)
						{
							initialkz[1]=i1;
							break;
						}
					}
					for(i1=initialkz[1];i1<N2;i1++)
					{
						if(total_xmin[1][i1]>x1)
						{
							if(i1>=1)
					    		finalkz[1]=i1-1;
							else
								finalkz[1]=i1;
							break;
						}
						else if(total_xmin[1][i1+1]==0)
						{
							finalkz[1]=i1;
							break;
						}
					}

                    
                    for(j1=0;j1<N1;j1++)
					{
		   				if(total_ymax[1][j1]>y0)
						{
							initialwz[1]=j1;
							break;
						}
					}
					for(j1=initialwz[1];j1<N1;j1++)
					{
						if(total_ymin[1][j1]>y1)
						{
							if(j1>=1)
					    		finalwz[1]=j1-1;
							else
								finalwz[1]=j1;
							break;
						}
						else if(total_ymin[1][j1+1]==0)
						{
							finalwz[1]=j1;
							break;
						}	
					}


/*以下确定第k-2行椭圆中与此椭圆定位有关的所有椭圆，第initialkz[0]列到第finalkz[0]列，第initialwz[0]行到第finalwz[0]行****/
					if(k>=2)
					{
						for(i1=0;i1<N2;i1++)
						{
	    	   				if(total_xmax[0][i1]>x0)
							{
		    					initialkz[0]=i1;
	    						break;
							}
						}
	    				for(i1=initialkz[0];i1<=N2;i1++)
						{
	    					if(total_xmin[0][i1]>x1)
							{
	    						if(i1>=1)
	    				    		finalkz[0]=i1-1;
   	    						else
			    					finalkz[0]=i1;
	 		    				break;
							}
		    				else if(total_xmin[0][i1+1]==0)
							{
    							finalkz[0]=i1;
	    						break;
							}
						}

                    
                        for(j1=0;j1<N1;j1++)
						{
    		   				if(total_ymax[0][j1]>y0)
							{
		    					initialwz[0]=j1;
		    					break;
							}
						}
    					for(j1=initialwz[0];j1<N1;j1++)
						{
	    					if(total_ymin[0][j1]>y1)
							{
		    					if(j1>=1)
		    			    		finalwz[0]=j1-1;
			    				else
		    						finalwz[0]=j1;
		    					break;
							}
		    				else if(total_ymin[0][j1+1]==0)
							{
		    					finalwz[0]=j1;
		    					break;
							}	
						}
					}

                    firstsettlez(elliptic,k,j,i);
					
					(*p0).zmin=(*p0).center_z-sqrt((*p0).C);
					(*p0).zmax=(*p0).center_z+sqrt((*p0).C);

			

				
				}
///*************************************z的第一次定位结束************************************************************///		    

///**************************************第二次定位***************************************************************///
				do
				{//dx,dy,dz为用此方法当然椭球在各个方向所能移动的距离
                    if(k>0)
	            	    dz=secondsettlez(elliptic,k,j,i);
			    	 else
			    		 dz=0;
			    	 if(j>0)
			    		 dy=secondsettley(elliptic,k,j,i);
		    		 else
		    			 dy=0;
                     if(i>0)
			    		 dx=secondsettlex(elliptic,k,j,i);
		    		 else
		    			 dx=0;
					
		    		 if(dz>=dy&&dz>=dx)
					 {
		    			(*p0).center_z=(*p0).center_z-dz;
		    			(*p0).zmin=(*p0).zmin-dz;
		    			(*p0).zmax=(*p0).zmax-dz;
					 }
		    		 else if(dy>=dz&&dy>=dx)
					 {
		    			 (*p0).center_y=(*p0).center_y-dy;
		    			(*p0).ymin=(*p0).ymin-dy;
			    		(*p0).ymax=(*p0).ymax-dy;
					 }
					 else if(dx>=dz&&dx>=dy)
					 {
		    			 (*p0).center_x=(*p0).center_x-dx;
		    			(*p0).xmin=(*p0).xmin-dx;
			    		(*p0).xmax=(*p0).xmax-dx;
					 }
				}while(dz>errort||dy>errort||dx>errort);
			
//printf("%d %d %d %lf %lf %lf\n",k,j,i,(*p0).center_x,(*p0).center_y,(*p0).center_z);

				 

				                                    
///*****************************************************二次定位结束*******************************************///

				 	{//当然椭球若满足插入的条件，则插入上一层或上一行
						double limitz;
						double limity;
						limitz=0;
						limity=0;

						if(k>0)
						{
							for(j1=initialwz[1];j1<=finalwz[1];j1++)
					    		for(i1=initialkz[1];i1<=finalkz[1];i1++)
						    		if(elliptic[1][j1][i1].zmax>limitz)
                                        limitz=elliptic[1][j1][i1].zmax;
		

				    		if((*p0).zmax<limitz+b_aver/5)
							{
					    		if(insertz(p0,elliptic[1]))
								{
				 	    	    	norm4=0;
								}
								else
									norm4=1;
							}
						}

						if(norm4==1)
						{
							if(j>0)
							{
					    		for(i1=initialy[1];i1<=finaly[1];i1++)
									if(elliptic[2][j-1][i1].ymax>limity)
										limity=elliptic[2][j-1][i1].ymax;

					    		if((*p0).ymax<limity+b_aver/5)
								{
									//mj++;
					    			if(inserty(p0,elliptic[2][j-1]))
									{
					    	    		norm5=0;
									}
									else
										norm5=1;
								}
							}
						}
								
					}
					//printf("norm4 norm5 %d,%d\n",norm4,norm5);
				{//当前椭球插入到链表
					(*p2).a=elliptic[2][j][i].a;
					(*p2).al=elliptic[2][j][i].al;
                    (*p2).b=elliptic[2][j][i].b;
					(*p2).bl=elliptic[2][j][i].bl;
					(*p2).c=elliptic[2][j][i].c;
					(*p2).cl=elliptic[2][j][i].cl;
					(*p2).center_x=elliptic[2][j][i].center_x;
					(*p2).center_y=elliptic[2][j][i].center_y;
					(*p2).center_z=elliptic[2][j][i].center_z;
					(*p2).alpha1=elliptic[2][j][i].alpha1;
					(*p2).alpha2=elliptic[2][j][i].alpha2;
					(*p2).alpha3=elliptic[2][j][i].alpha3;
					(*p2).bata1=elliptic[2][j][i].bata1;
					(*p2).bata2=elliptic[2][j][i].bata2;
					(*p2).bata3=elliptic[2][j][i].bata3;
					(*p2).gama1=elliptic[2][j][i].gama1;
					(*p2).gama2=elliptic[2][j][i].gama2;
					(*p2).gama3=elliptic[2][j][i].gama3;
					(*p2).A=elliptic[2][j][i].A;
					(*p2).B=elliptic[2][j][i].B;
					(*p2).C=elliptic[2][j][i].C;

	//	printf("%d %d %d %d %d %lf %lf %lf\n",norm4,norm5,k,j,i,(*p0).center_z,(*p0).center_y,(*p0).center_x);

					if(norm4==0||norm5==0)
					{//若当前椭球插入上一层或者上一行则重新生成当然椭球
						i=i-1;
						continue;
					}
				}

				//当当前椭球的xmax>=rightbar1，则开始生成下一行；当前行的所有椭球的ymax>=formbar1,则开始生成下一层；
				//若当前层的所有椭球的zmax〉=upbar1,则生成结束
				if(elliptic[2][j][i].zmax<upbar1)
					norm1=0;
				if(elliptic[2][j][i].ymax<formbar1)
            	   norm2=0;
            	if(elliptic[2][j][i].xmax>=rightbar1)
				{
            	   break;
				}
				
			}
			if(norm2==1)
			{//printf("%lf\n",total_ymin[0][0]);
				
				break;
			}
		}

	
			
		    for(j1=0;j1<N1;j1++)
	    	    for(i1=0;i1<N2;i1++)
				{
	        		replace(&elliptic[0][j1][i1],&elliptic[1][j1][i1]);
	        		replace(&elliptic[1][j1][i1],&elliptic[2][j1][i1]);
				}

				//每一排中，每一行以及每一列，x,y总的最小，最大值
			for(k1=0;k1<=1;k1++)
			{
				for(j1=0;j1<N1;j1++)
				{
					for(i1=0;i1<N2;i1++)
					{
						if(elliptic[k1][j1][i1].al!=0)
						{
							total_ymin[k1][j1]=elliptic[k1][j1][i1].ymin;
							total_ymax[k1][j1]=elliptic[k1][j1][i1].ymax;
							break;
						}
					}
			
     				for(;i1<N2;i1++)
					{
						if(elliptic[k1][j1][i1].ymin<total_ymin[k1][j1]&&elliptic[k1][j1][i1].al!=0)
							total_ymin[k1][j1]=elliptic[k1][j1][i1].ymin;
						if(elliptic[k1][j1][i1].ymax>total_ymax[k1][j1]&&elliptic[k1][j1][i1].al!=0)
                            total_ymax[k1][j1]=elliptic[k1][j1][i1].ymax;
					}
			
				}
			}

			for(k1=0;k1<=1;k1++)
				for(i1=0;i1<N2;i1++)
				{
					for(j1=0;j1<N1;j1++)
					{
						if(elliptic[k1][j1][i1].al!=0)
						{
							total_xmin[k1][i1]=elliptic[k1][j1][i1].xmin;
							total_xmax[k1][i1]=elliptic[k1][j1][i1].xmax;
					    	break;
						}
					}
					for(;j1<N1;j1++)
					{
						if(elliptic[k1][j1][i1].xmin<total_xmin[k1][i1]&&elliptic[k1][j1][i1].al!=0)
							total_xmin[k1][i1]=elliptic[k1][j1][i1].xmin;
						if(elliptic[k1][j1][i1].xmax>total_xmax[k1][i1]&&elliptic[k1][j1][i1].al!=0)
							total_xmax[k1][i1]=elliptic[k1][j1][i1].xmax;

					}
				}


			if(norm1==1)
				break;
	}

	(*p2).next=NULL;

	for(k1=0;k1<=1;k1++)
	{
		free(total_xmin[k1]);
		free(total_xmax[k1]);
		free(total_ymin[k1]);
		free(total_ymax[k1]);
	}
	for(k1=0; k1<=2;k1++)
	{
		free(elliptic[k1]);
	}

	return 1;

//printf("f e c %d %d %d\n",f_num,e_num,c_num);
	
};	

int EllipseGen::uniform(ifstream &input_stream, string outfile)
{
	int compactnumber=0,needednumber=0;
	double compactrate=0,rate=0,ratemin,ratemax;
	int i=0,k=0;    
	double p1,p2,q1,q2,r1,r2,p11,p22,q11,q22,r11,r22;
	struct parameter *p0;

//均匀分布体积百分含量的读入（<紧凑分布体积百分含量）
	char rline[250];
	istrstream istr(get_line(input_stream,rline));
    istr >> neededrate;

    p0=head;
    do
	{	
		if((*p0).center_x>=leftbar&&(*p0).center_x<=rightbar
		   	    	&&(*p0).center_y>=backbar&&(*p0).center_y<=formbar
			           	&&(*p0).center_z>=downbar&&(*p0).center_z<=upbar)
		{
			compactnumber++;
			compactrate=compactrate+volume((*p0).center_x,(*p0).center_y,(*p0).center_z,p0);
		}

		p0=(*p0).next;

	}while(p0!=NULL);

	compactrate=compactrate/((rightbar-leftbar)*(upbar-downbar)*(formbar-backbar));
	ratemax=compactrate;
	ratemin=0;

    //紧凑分布体积百分含量
   printf("紧凑分布体积百分含量 %lf\n",compactrate);
   
    if(neededrate<compactrate)
		printf("The needed rate can be reached.\n");
    else 
    {
    	 printf("The needed rate is higher than that can be reached.");
    	 return 0;
    }


	p1=(leftbar+rightbar)/2-0.5*pow((neededrate),0.34)*(rightbar-leftbar);
	p11=leftbar;
	p2=p1+pow((neededrate),0.34)*(rightbar-leftbar);
	p22=rightbar;

	q1=(formbar+backbar)/2-0.5*pow((neededrate),0.34)*(formbar-backbar);;
	q11=backbar;
	q2=q1+pow((neededrate),0.34)*(formbar-backbar);
	q22=formbar;

	r1=(downbar+upbar)/2-0.5*pow((neededrate),0.34)*(upbar-downbar);
	r11=downbar;
	r2=r1+pow((neededrate),0.34)*(upbar-downbar);
	r22=upbar;



	do
	{
		rate=0;
		leftbar1=(p1+p11)/2;
		rightbar1=(p2+p22)/2;
		backbar1=(q1+q11)/2;
		formbar1=(q2+q22)/2;
		downbar1=(r1+r11)/2;
		upbar1=(r2+r22)/2;

        p0=head;
    	do
		{	
	    	if((*p0).center_x>leftbar1&&(*p0).center_x<rightbar1
		   	    	&&(*p0).center_y>backbar1&&(*p0).center_y<formbar1
			           	&&(*p0).center_z>downbar1&&(*p0).center_z<upbar1)
			{
				double x1,y1,z1;
				x1=((*p0).center_x-leftbar1)*(rightbar-leftbar)/(rightbar1-leftbar1)+leftbar;
		    	y1=((*p0).center_y-backbar1)*(formbar-backbar)/(formbar1-backbar1)+backbar;
		    	z1=((*p0).center_z-downbar1)*(upbar-downbar)/(upbar1-downbar1)+downbar;
		    	rate=rate+volume(x1,y1,z1,p0);//////
		
			}
        
	    	p0=(*p0).next;

		}while(p0!=NULL);

	    rate=rate/((rightbar-leftbar)*(formbar-backbar)*(upbar-downbar));

		if(rate<neededrate)
		{
			p1=leftbar1;
			p2=rightbar1;
			q1=backbar1;
			q2=formbar1;
			r1=downbar1;
			r2=upbar1;
			ratemin=rate;
		}
		else
		{
			p11=leftbar1;
			p22=rightbar1;
			q11=backbar1;
			q22=formbar1;
			r11=downbar1;
			r22=upbar1;
			ratemax=rate;

		}

		double V_rate;
		V_rate=(p22*r22*q22-p2*q2*r2)/((rightbar-leftbar)*(formbar-backbar)*(upbar-downbar));//二分法中相邻俩个长方体区域体积比之差
		if(V_rate<0.00000000001)
			break;
	//	printf("V_rate r rmin rmax %lf %lf %lf %lf\n",V_rate,rate,ratemin,ratemax);
		if((ratemax-neededrate)*(ratemin-neededrate)>0)/////
			break;

		
	}while(fabs(rate-neededrate)>0.0001);

	
    printf("满足要求的椭球的比率 rate %lf\n",rate);
	realrate=rate;

	needednumber=0;
    p0=head;
   	do
	{	
	   	if((*p0).center_x>=leftbar1&&(*p0).center_x<=rightbar1
		   	   	&&(*p0).center_y>=backbar1&&(*p0).center_y<=formbar1
		           	&&(*p0).center_z>=downbar1&&(*p0).center_z<=upbar1)
		{
		   	needednumber++;
		}
        
	   	p0=(*p0).next;

	}while(p0!=NULL);

	/*   printf the bata. */

	ofstream output_stream;
	output_stream.open(outfile.c_str(),ios::out);

    if(!output_stream)
	{
		cout <<"Cannot open file.\n" <<outfile <<endl;
		return 0;
	}

	output_stream <<"%单胞内颗粒个数和体积分数" <<endl;
	output_stream << needednumber << " " << rate << endl;

    //输出单胞位置和尺寸
	output_stream <<"%单胞起始点坐标" <<endl;
	output_stream <<leftbar << " " <<backbar << " " <<downbar <<endl;
	output_stream <<"%单胞长宽高" <<endl;
	output_stream <<rightbar << " " <<formbar << " " <<upbar <<endl;

	//输出椭球15个参数（参数方程：椭球中心坐标（3个）、尺寸（3个）、方位角（9个））
	output_stream <<"%输出椭球15个参数:椭球中心坐标、尺寸、方位角" <<endl;


	do
	{
		if((*p0).center_x>=leftbar1&&(*p0).center_x<=rightbar1
			&&(*p0).center_y>=backbar1&&(*p0).center_y<=formbar1
		   &&(*p0).center_z>=downbar1&&(*p0).center_z<=upbar1)
		{
			double X,Y,Z;
			X=((*p0).center_x-leftbar1)*(rightbar-leftbar)/(rightbar1-leftbar1)+leftbar;
			Y=((*p0).center_y-backbar1)*(formbar-backbar)/(formbar1-backbar1)+backbar;
			Z=((*p0).center_z-downbar1)*(upbar-downbar)/(upbar1-downbar1)+downbar;
            output_stream  <<X <<" " <<Y <<" " <<Z <<" "
				<<(*p0).a <<" " <<(*p0).b <<" " <<(*p0).c <<" "
				<<(*p0).alpha1 <<" " <<(*p0).alpha2 <<" " <<(*p0).alpha3 <<" "
                <<(*p0).bata1 <<" " <<(*p0).bata2 <<" " <<(*p0).bata3 <<" "
				<<(*p0).gama1 <<" " <<(*p0).gama2 <<" " <<(*p0).gama3 <<endl;
		}

		p0=(*p0).next;

	}while(p0!=NULL);

	output_stream.close();
	free(head);

	return 1;
};




int EllipseGen::no_uniform(int distribution_sign2,ifstream& input_stream,string out_file, string data_file)
{
	struct parameter *p;
	struct basicparameter *ellipsoid,*p0;
	int i,j,k=0,k1=0;
	double sitax,miux,sitay,miuy,sitaz,miuz;
	double arx,ary,arz;
	double leftbar1,rightbar1,backbar1,formbar1,downbar1,upbar1;
    double rate=0;
    int compactnumber=0,numbermax=0,needednumber=0,uninumber=0;
    int mx=0,my=0,mz=0,m=0;
    double *sequencex,*sequencey,*sequencez;
    double gamax,gamay,gamaz;
    double *ratiox,*ratioy,*ratioz;
	vector<int> orderx,ordery,orderz,order;
	char rline[250];

	if(distribution_sign2==0)
   {
     	gamax=(rightbar-leftbar);
    	gamay=(formbar-backbar);
		gamaz=(upbar-downbar);
	}

	if(distribution_sign2==1)
   {
		istrstream istr2(get_line(input_stream,rline));
		istr2 >>sitax >>miux >>sitay >>miuy >>sitaz >>miuz ;
     	gamax=sqrt(2*pi)*sitax;
    	gamay=sqrt(2*pi)*sitax;
		gamaz=sqrt(2*pi)*sitax;
	}

	if(distribution_sign2==2)
	{
		istrstream istr2(get_line(input_stream,rline));
		istr2 >> arx >>ary >>arz ;
		gamax=1/arx;
		gamay=1/ary;
		gamaz=1/arz;

	}

	//非均匀分布体积百分含量的读入（<紧凑分布体积百分含量）
	
	istrstream istr1(get_line(input_stream,rline));
    istr1 >> neededrate ;

	p=head;
	do
	 {	
		 /* if((*p).center_x-sqrt((*p).A)>leftbar&&(*p).center_x+sqrt((*p).A)<rightbar
		  	&&(*p).center_y-sqrt((*p).B)>backbar&&(*p).center_y+sqrt((*p).B)<formbar
		   	&&(*p).center_z-sqrt((*p).C)>downbar&&(*p).center_z+sqrt((*p).C)<upbar)*/
	   	if((*p).center_x>=leftbar&&(*p).center_x<=rightbar
		  	&&(*p).center_y>=backbar&&(*p).center_y<=formbar
		   	&&(*p).center_z>=downbar&&(*p).center_z<=upbar)
		{
	   		compactnumber++;
		}
   		p=(*p).next;
	 }while(p!=NULL);

	 ellipsoid=(struct basicparameter *)malloc(compactnumber*sizeof(struct basicparameter));
     p=head;
	 i=0;
	do
	 {	
		 /*if((*p).center_x-sqrt((*p).A)>leftbar&&(*p).center_x+sqrt((*p).A)<rightbar
		  	&&(*p).center_y-sqrt((*p).B)>backbar&&(*p).center_y+sqrt((*p).B)<formbar
		   	&&(*p).center_z-sqrt((*p).C)>downbar&&(*p).center_z+sqrt((*p).C)<upbar)*/
	   	if((*p).center_x>=leftbar&&(*p).center_x<=rightbar
		  	&&(*p).center_y>=backbar&&(*p).center_y<=formbar
		   	&&(*p).center_z>=downbar&&(*p).center_z<=upbar)
		{
	   		ellipsoid[i].a=(*p).a;
			ellipsoid[i].b=(*p).b;
			ellipsoid[i].c=(*p).c;
			ellipsoid[i].al=(*p).al;
			ellipsoid[i].bl=(*p).bl;
			ellipsoid[i].cl=(*p).cl;
			ellipsoid[i].alpha1=(*p).alpha1;
			ellipsoid[i].alpha2=(*p).alpha2;
			ellipsoid[i].alpha3=(*p).alpha3;
			ellipsoid[i].bata1=(*p).bata1;
			ellipsoid[i].bata2=(*p).bata2;
			ellipsoid[i].bata3=(*p).bata3;
			ellipsoid[i].gama1=(*p).gama1;
			ellipsoid[i].gama2=(*p).gama2;
			ellipsoid[i].gama3=(*p).gama3;
			ellipsoid[i].center_x=(*p).center_x;
			ellipsoid[i].center_y=(*p).center_y;
			ellipsoid[i].center_z=(*p).center_z;
			ellipsoid[i].A=(*p).A;
			ellipsoid[i].B=(*p).B;
			ellipsoid[i].C=(*p).C;
			i++;
		}
   		p=(*p).next;
	 }while(p!=NULL);

    free(head);

	leftbar1=leftbar;
	rightbar1=rightbar;
	backbar1=backbar;
	formbar1=formbar;
	downbar1=downbar;
	upbar1=upbar;

		double *center_x,*center_y,*center_z;
		k1++;
		rate=0;
		uninumber=0;

		 
        for(i=0;i<compactnumber;i++)
		 {
			p0=&ellipsoid[i];
	    	if((*p0).center_x>=leftbar1&&(*p0).center_x<=rightbar1
		    	&&(*p0).center_y>=backbar1&&(*p0).center_y<=formbar1
		    	&&(*p0).center_z>=downbar1&&(*p0).center_z<=upbar1)

			{
	    		uninumber++;
			}
		}

	//------------------------------------------------------------------------------------------------------------------------
	//生成的紧凑椭球结果输出(080616韩非增加)
	//ofstream out_compact;
	//out_compact.open("CompactResult.dat");
 //   if(!out_compact)
	//{
	//	cout <<"Cannot open file.\n" <<"CompactResult.dat" <<endl;
	//	return 0;
	//}
	//out_compact <<"%单胞内颗粒个数和体积分数" <<endl;
	//out_compact << compactnumber << endl;
 //   //输出单胞位置和尺寸
	//out_compact <<"%单胞起始点坐标" <<endl;
	//out_compact <<leftbar << " " <<backbar << " " <<downbar <<endl;
	//out_compact <<"%单胞长宽高" <<endl;
	//out_compact <<rightbar << " " <<formbar << " " <<upbar <<endl;
	////输出椭球15个参数（参数方程：椭球中心坐标（3个）、尺寸（3个）、方位角（9个））
	//out_compact <<"%输出椭球15个参数:椭球中心坐标、尺寸、方位角" <<endl;
	//for(j=0;j<compactnumber;j++)
	//{
	//	double X,Y,Z;
	//	p0=&ellipsoid[j];
	//	X=((*p0).center_x-leftbar1)*(rightbar-leftbar)/(rightbar1-leftbar1)+leftbar;
	//	Y=((*p0).center_y-backbar1)*(formbar-backbar)/(formbar1-backbar1)+backbar;
	//	Z=((*p0).center_z-downbar1)*(upbar-downbar)/(upbar1-downbar1)+downbar;
 //       out_compact  <<X <<" " <<Y <<" " <<Z <<" "
	//			<<(*p0).a <<" " <<(*p0).b <<" " <<(*p0).c <<" "
	//			<<(*p0).alpha1 <<" " <<(*p0).alpha2 <<" " <<(*p0).alpha3 <<" "
 //               <<(*p0).bata1 <<" " <<(*p0).bata2 <<" " <<(*p0).bata3 <<" "
	//			<<(*p0).gama1 <<" " <<(*p0).gama2 <<" " <<(*p0).gama3 <<endl;
	//}
	//out_compact.close();
	//------------------------------------------------------------------------------------------------------------------------

    	center_x=(double *)malloc(uninumber*sizeof(double));
    	center_y=(double *)malloc(uninumber*sizeof(double));
    	center_z=(double *)malloc(uninumber*sizeof(double));

        for(i=0;i<uninumber;i++)
		{
			p0=&ellipsoid[i];
    		if((*p0).center_x>=leftbar1&&(*p0).center_x<=rightbar1
	    		&&(*p0).center_y>=backbar1&&(*p0).center_y<=formbar1
		    	&&(*p0).center_z>=downbar1&&(*p0).center_z<=upbar1)

			{
	            center_x[i]=(*p0).center_x ; 	
				center_y[i]=(*p0).center_y;
				center_z[i]=(*p0).center_z;
			}	
		 }

    	sequencex=(double *)malloc(uninumber*sizeof(double));
        sequencey=(double *)malloc(uninumber*sizeof(double));
    	sequencez=(double *)malloc(uninumber*sizeof(double));
        for(i=0;i<uninumber;i++)
		{
       	  sequencex[i]=unifrnd(0.0,1.0,&seqx);
		}
        for(i=0;i<uninumber;i++)
		{
    	  sequencey[i]=unifrnd(0.0,1.0,&seqy);
		}
    	for(i=0;i<uninumber;i++)
		{
    		sequencez[i]=unifrnd(0.0,1.0,&seqz);
		}


    

  
        ratiox=(double *)malloc(uninumber*sizeof(double));
        ratioy=(double *)malloc(uninumber*sizeof(double));
    	ratioz=(double *)malloc(uninumber*sizeof(double));

		if(distribution_sign2==0)
		{
        	for(i=0;i<uninumber;i++)
			{
         	  ratiox[i]=sequencex[i]*(rightbar-leftbar);
        	  ratioy[i]=sequencey[i]*(formbar-backbar);
              ratioz[i]=sequencez[i]*(upbar-downbar);
		  
			}
		}

        if(distribution_sign2==1)
		{
        	for(i=0;i<uninumber;i++)
			{
         	  ratiox[i]=sequencex[i]/normfunction(miux,sitax,center_x[i]);
        	  ratioy[i]=sequencey[i]/normfunction(miuy,sitay,center_y[i]);
              ratioz[i]=sequencez[i]/normfunction(miuz,sitaz,center_z[i]);
		  
			}
		}
 

       if(distribution_sign2==2)
	   {
        	for(i=0;i<uninumber;i++)
			{
        	  ratiox[i]=sequencex[i]/expfunction(arx,center_x[i]);
        	  ratioy[i]=sequencey[i]/expfunction(ary,center_y[i]);
    		  ratioz[i]=sequencez[i]/expfunction(arz,center_z[i]);
	      
			} 
	   }
  
        free(sequencex);
        free(sequencey);
    	free(sequencez);
    	
   for(i=0;i<uninumber;i++)
   {
     orderx.push_back(i);
	 ordery.push_back(i);
	 orderz.push_back(i);
   }
		//对比值进行排序
  for(i=0;i<uninumber;i++)
  {   
	  double cf;
	  int cd;

	  for(j=i+1;j<uninumber;j++)
	  {
	    if(ratiox[i]>ratiox[j])
		{
		  cf=ratiox[j];
		  ratiox[j]=ratiox[i];
		  ratiox[i]=cf;
		  cd=orderx[j];
		  orderx[j]=orderx[i];
		  orderx[i]=cd;
		}
		if(ratioy[i]>ratioy[j])
		{
		  cf=ratioy[j];
		  ratioy[j]=ratioy[i];
		  ratioy[i]=cf;
		  cd=ordery[j];
		  ordery[j]=ordery[i];
		  ordery[i]=cd;
		}
        	if(ratioz[i]>ratioz[j])
		{
		  cf=ratioz[j];
		  ratioz[j]=ratioz[i];
		  ratioz[i]=cf;
		  cd=orderz[j];
		  orderz[j]=orderz[i];
		  orderz[i]=cd;
		}
	  }
  }

//满足αifi(ξi(1))≥η 的椭球序号保留
       vector<int>::iterator iter=orderx.begin();
        for(i=0;i<uninumber;i++)
		{
    	  if(ratiox[i]<gamax)
		  {
    		 continue;
		  }
		  else
		  {
			  iter+=i;
			  orderx.erase(iter,orderx.end());
			  break;
		  }
		}

		iter=ordery.begin();
        for(i=0;i<uninumber;i++)
		{
    	  if(ratioy[i]<gamay)
    	 {
    		 continue;
		  }
		  else
		 {
			  iter+=i;
			  ordery.erase(iter,ordery.end());
			  break;
		  }
		}

		iter=orderz.begin();
        for(i=0;i<uninumber;i++)
		{
    	   if(ratioz[i]<gamaz)
    	 {
    		 continue;
		  }
		   else
		 {
			  iter+=i;
			  orderz.erase(iter,orderz.end());
			  break;
		  }
		}
	
//二分法循环确定满足要求的椭球
	mx=(int)orderx.size();
	my=(int)ordery.size();
	mz=(int)orderz.size();

	k1=0;                        //循环的步数
	int mzmin,mzmax;
	mzmin=0;
	mzmax=mz;
     do
	 {
		 rate=0;
		 k1++;

		 if( (int)order.size()>0 )
			 order.erase( order.begin(),order.end() );

        for(k=0;k<mz;k++)
		{
    	  for(j=0;j<my;j++)
		  {
    		  if(ordery[j]==orderz[k])
			  {
        		for(i=0;i<mx;i++)
				{
        		  if(orderx[i]==ordery[j])
				  {
        			order.push_back(orderx[i]);	
					
        			rate=rate
						+volume(center_x[orderx[i]],center_y[orderx[i]],center_z[orderx[i]],&(ellipsoid[orderx[i]]))
						/((rightbar-leftbar)*(formbar-backbar)*(upbar-downbar));
					
        			goto label1;
				  }
				}
	    		goto label1;
			  }
		  }
label1: ;
		}
        needednumber=(int)order.size();
        
        if(k1==1)
		{
			hout << "最大比率：" << rate << endl;
    		if(neededrate<=rate)
               printf("The needed rate can be reached.\n");
            else 
			{
        		printf("The needed rate is higher than that can be reached.");
        		return 0;
			}
		}

		if(rate<neededrate)
		{
			mzmin=mz;
		}
		else
		{
			mzmax=mz;
		}

		if( mzmax-mzmin<=1 )
			break;
		else
    		mz=(int)( (mzmin+mzmax)/2.0 );

	}while(1);

//	 printf("实际比率rate=%lf\n",rate);

	//输入到样本数据文件
	ofstream out(data_file.c_str(),ios::app);
	out <<"%单胞内颗粒个数和体积分数" << endl;
	out << needednumber << " " << rate <<endl;
	out.close();

	realrate=rate;

 	free(center_x);
   	free(center_y);
   	free(center_z);
	free(ratiox);
    free(ratioy);
    free(ratioz);

	ofstream output_stream;
	output_stream.open(out_file.c_str(),ios::out);

    if(!output_stream)
	{
		cout <<"Cannot open file.\n" <<out_file <<endl;
		return 0;
	}

	output_stream <<"%单胞内颗粒个数和体积分数" <<endl;
	output_stream << needednumber << " " << rate << endl;

    //输出单胞位置和尺寸
	output_stream <<"%单胞起始点坐标" <<endl;
	output_stream <<leftbar << " " <<backbar << " " <<downbar <<endl;
	output_stream <<"%单胞长宽高" <<endl;
	output_stream <<rightbar << " " <<formbar << " " <<upbar <<endl;

 //输出椭球15个参数（参数方程：椭球中心坐标（3个）、尺寸（3个）、方位角（9个））
	output_stream <<"%输出椭球15个参数:椭球中心坐标、尺寸、方位角" <<endl;

	for(j=0;j<needednumber;j++)
	{
		double X,Y,Z;
		p0=&ellipsoid[order[j]];
		X=((*p0).center_x-leftbar1)*(rightbar-leftbar)/(rightbar1-leftbar1)+leftbar;
		Y=((*p0).center_y-backbar1)*(formbar-backbar)/(formbar1-backbar1)+backbar;
		Z=((*p0).center_z-downbar1)*(upbar-downbar)/(upbar1-downbar1)+downbar;
        output_stream  <<X <<" " <<Y <<" " <<Z <<" "
				<<(*p0).a <<" " <<(*p0).b <<" " <<(*p0).c <<" "
				<<(*p0).alpha1 <<" " <<(*p0).alpha2 <<" " <<(*p0).alpha3 <<" "
                <<(*p0).bata1 <<" " <<(*p0).bata2 <<" " <<(*p0).bata3 <<" "
				<<(*p0).gama1 <<" " <<(*p0).gama2 <<" " <<(*p0).gama3 <<endl;
	}

	output_stream.close();

	return 1;
}





void EllipseGen::firstsettley(struct parameter0 **p,int j,int i)
{
	double ymax=0;
	int i1;
	int norm=0;

	for(i1=initialy[1];i1<=finaly[1];i1++)
		if(p[j-1][i1].ymax>ymax)
		{
			norm=1;
			ymax=p[j-1][i1].ymax;
		}
	if(j>=2)
	{
		for(i1=initialy[0];i1<=finaly[0];i1++)
     		if(p[j-2][i1].ymax>ymax)
			{
				norm=1;
     			ymax=p[j-2][i1].ymax;
			}
	}

	if(norm==0)
		ymax=p[j-1][0].ymax;
	p[j][i].center_y=ymax+sqrt(p[j][i].B);
};



void EllipseGen::firstsettlez(struct parameter0 ***p,int k,int j,int i)
{
	double zmax=0;
	int i1,j1;
	int norm=0;

	for(j1=initialwz[1];j1<=finalwz[1];j1++)
		for(i1=initialkz[1];i1<=finalkz[1];i1++)
		{
			norm=1;
			if(p[1][j1][i1].zmax>zmax)
				zmax=p[1][j1][i1].zmax;
		}
	if(k>=2)
	{
    	for(j1=initialwz[0];j1<=finalwz[0];j1++)
    		for(i1=initialkz[0];i1<=finalkz[0];i1++)
			{
				norm=1;
    			if(p[0][j1][i1].zmax>zmax)
    				zmax=p[0][j1][i1].zmax;
			}
	}

	if(norm==0)
		zmax=p[1][0][0].zmax;
	p[2][j][i].center_z=zmax+sqrt(p[2][j][i].C);
};

 /**********************************************************************/
  /*                    z的第二次定位                                */
 /**********************************************************************/

double EllipseGen::secondsettlez(struct parameter0 ***p0,int k,int j,int i)
{
	double r1,r;
	int i1,j1;
	double xmin,xmax,ymin,ymax;
	int norm=0;
	

	r=p0[2][j][i].zmin-p0[1][0][0].zmin;

	for(j1=initialwz[1];j1<=finalwz[1];j1++)
		for(i1=initialkz[1];i1<=finalkz[1];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[1][j1][i1].xmin);
			xmax=min(p0[2][j][i].xmax,p0[1][j1][i1].xmax);
			ymin=max(p0[2][j][i].ymin,p0[1][j1][i1].ymin);
			ymax=min(p0[2][j][i].ymax,p0[1][j1][i1].ymax);
			
			if(xmin<xmax&&ymin<ymax)
			{
				norm=1;
				r1=min(dminz(&p0[1][j1][i1],&p0[2][j][i]),p0[2][j][i].zmin-p0[1][j1][i1].zmin);
				if(r1<r)
					r=r1;	
			}
		}

	if(k>1)
	{
		for(j1=initialwz[0];j1<=finalwz[0];j1++)
		for(i1=initialkz[0];i1<=finalkz[0];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[0][j1][i1].xmin);
			xmax=min(p0[2][j][i].xmax,p0[0][j1][i1].xmax);
			ymin=max(p0[2][j][i].ymin,p0[0][j1][i1].ymin);
			ymax=min(p0[2][j][i].ymax,p0[0][j1][i1].ymax);
			if(xmin<xmax&&ymin<ymax)
			{
			    norm=1;
				r1=dminz(&p0[0][j1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}


	if(j>=1)
	{
		for(i1=initialy[1];i1<=finaly[1];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j-1][i1].xmin);
			xmax=min(p0[2][j][i].xmax,p0[2][j-1][i1].xmax);
			ymin=max(p0[2][j][i].ymin,p0[2][j-1][i1].ymin);
			ymax=min(p0[2][j][i].ymax,p0[2][j-1][i1].ymax);
			if(xmin<xmax&&ymin<ymax)
			{
			    norm=1;
				r1=dminz(&p0[2][j-1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(j>=2)
	{
		for(i1=initialy[0];i1<=finaly[0];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j-2][i1].xmin);
			xmax=min(p0[2][j][i].xmax,p0[2][j-2][i1].xmax);
			ymin=max(p0[2][j][i].ymin,p0[2][j-2][i1].ymin);
			ymax=min(p0[2][j][i].ymax,p0[2][j-2][i1].ymax);
			if(xmin<xmax&&ymin<ymax)
			{
			    norm=1;
				r1=dminz(&p0[2][j-2][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(i>=1)
	{
		for(i1=initialx;i1<=finalx;i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j][i1].xmin);
			xmax=min(p0[2][j][i].xmax,p0[2][j][i1].xmax);
			ymin=max(p0[2][j][i].ymin,p0[2][j][i1].ymin);
			ymax=min(p0[2][j][i].ymax,p0[2][j][i1].ymax);
			if(xmin<xmax&&ymin<ymax)
			{
			    norm=1;
				r1=dminz(&p0[2][j][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}
//printf("norm zr %d %lf\n",norm,r);
	
	if(norm==0)
		r=0;
	return(r);
};


double EllipseGen::dminz(struct parameter0 *p1,struct parameter0 *p2)
{
	struct coordinate point[21][21];
	double xmin,xmax,ymin,ymax;
	double r=0,r1,r2;
	int i,j,k1,k2;
	int norm1=0,norm2=0;

	xmin=max((*p1).xmin,(*p2).xmin)+error;
	xmax=min((*p1).xmax,(*p2).xmax)-error;
	ymin=max((*p1).ymin,(*p2).ymin)+error;
	ymax=min((*p1).ymax,(*p2).ymax)-error;
	r1=(*p2).zmin-(*p1).zmin;
	 r2=r1;
    
	for(j=0;j<21;j++)
		for(i=0;i<21;i++)
		{
			point[j][i].x=xmin+i*(xmax-xmin)/20;
			point[j][i].y=ymin+j*(ymax-ymin)/20;
		}

	for(j=0;j<21;j++)
	{
		for(i=0;i<21;i++)
		{
			if(judgez(p1,p2,point[j][i].x,point[j][i].y))
			{
	    		r=disz(p1,p2,point[j][i].x,point[j][i].y);
					
				if(r<r1)
				{
		    	   r1=r;
		    	   k1=j;
		    	   k2=i;
				}
			}
			else
				continue;
					
		}
	}


	return(r1);

};

double EllipseGen::disz(struct parameter0 *p1,struct parameter0 *p2,double x,double y)
{
	double r;
	double A1,B1,C1,D1,E1,F1;
	double A2,B2,C2,D2,E2,F2;
	A1=(pow((*p1).coe[5],2)-(*p1).coe[2]*(*p1).coe[0])/pow((*p1).coe[2],2);
	A2=(pow((*p2).coe[5],2)-(*p2).coe[2]*(*p2).coe[0])/pow((*p2).coe[2],2);
	B1=(pow((*p1).coe[4],2)-(*p1).coe[2]*(*p1).coe[1])/pow((*p1).coe[2],2);
    B2=(pow((*p2).coe[4],2)-(*p2).coe[2]*(*p2).coe[1])/pow((*p2).coe[2],2);
	C1=((*p1).coe[4]*(*p1).coe[5]-(*p1).coe[2]*(*p1).coe[3])/pow((*p1).coe[2],2);
	C2=((*p2).coe[4]*(*p2).coe[5]-(*p2).coe[2]*(*p2).coe[3])/pow((*p2).coe[2],2);
	D1=(*p1).coe[6]/(*p1).coe[2];
	D2=(*p2).coe[6]/(*p2).coe[2];
	E1=-(*p1).coe[5]/(*p1).coe[2];
	E2=-(*p2).coe[5]/(*p2).coe[2];
    F1=-(*p1).coe[4]/(*p1).coe[2];
	F2=-(*p2).coe[4]/(*p2).coe[2];

	r=-sqrt(A2*pow(x-(*p2).center_x,2)+B2*pow(y-(*p2).center_y,2)+2*C2*(x-(*p2).center_x)*(y-(*p2).center_y)+D2)
		        +E2*(x-(*p2).center_x)+F2*(y-(*p2).center_y)+(*p2).center_z
	   -sqrt(A1*pow(x-(*p1).center_x,2)+B1*pow(y-(*p1).center_y,2)+2*C1*(x-(*p1).center_x)*(y-(*p1).center_y)+D1)
	            -E1*(x-(*p1).center_x)-F1*(y-(*p1).center_y)-(*p1).center_z;
	return(r);
}

int EllipseGen::judgez(struct parameter0 *p1,struct parameter0 *p2,double x,double y)
{
	double A1,B1,C1,D1;
	double A2,B2,C2,D2;
	int norm;

	A1=(*p1).coe[2]*(*p1).coe[0]-pow((*p1).coe[5],2);
	B1=(*p1).coe[1]*(*p1).coe[2]-pow((*p1).coe[4],2);
	C1=(*p1).coe[2]*(*p1).coe[3]-(*p1).coe[4]*(*p1).coe[5];
	D1=(*p1).coe[6]*(*p1).coe[2];

	A2=(*p2).coe[2]*(*p2).coe[0]-pow((*p2).coe[5],2);
	B2=(*p2).coe[1]*(*p2).coe[2]-pow((*p2).coe[4],2);
	C2=(*p2).coe[2]*(*p2).coe[3]-(*p2).coe[4]*(*p2).coe[5];
	D2=(*p2).coe[6]*(*p2).coe[2];

	if(A1*pow(x-(*p1).center_x,2)+B1*pow(y-(*p1).center_y,2)+C1*2*(x-(*p1).center_x)*(y-(*p1).center_y)<=D1
	    &&A2*pow(x-(*p2).center_x,2)+B2*pow(y-(*p2).center_y,2)+C2*2*(x-(*p2).center_x)*(y-(*p2).center_y)<=D2)
		norm=1;
	else
		norm=0;

	return(norm);
}


/*******************************************************************************************************/
/*                              z定位完毕                             */
 /******************************************************************************************************/




 /**********************************************************************/
  /*                    y的第二次定位                                */
 /**********************************************************************/

double EllipseGen::secondsettley(struct parameter0 ***p0,int k,int j,int i)
{
	double r1,r;
	int i1,j1;
	double xmin,xmax,zmin,zmax;
	int norm=0;


	if(j>=1)
    	r=p0[2][j][i].ymin-p0[2][j-1][0].ymin;
	else
		r=p0[2][j][i].ymin;

	

	for(j1=initialwz[1];j1<=finalwz[1];j1++)
		for(i1=initialkz[1];i1<=finalkz[1];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[1][j1][i1].xmin)+error;
			xmax=min(p0[2][j][i].xmax,p0[1][j1][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[1][j1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[1][j1][i1].zmax)-error;
			if(xmin<xmax&&zmin<zmax)
			{
				norm=1;
				r1=dminy(&p0[1][j1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}

	if(k>1)
	{
		for(j1=initialwz[0];j1<=finalwz[0];j1++)
		for(i1=initialkz[0];i1<=finalkz[0];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[0][j1][i1].xmin)+error;
			xmax=min(p0[2][j][i].xmax,p0[0][j1][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[0][j1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[0][j1][i1].zmax)-error;
			if(xmin<xmax&&zmin<zmax)
			{
			    norm=1;
				r1=dminy(&p0[0][j1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}


	if(j>=1)
	{
		for(i1=initialy[1];i1<=finaly[1];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j-1][i1].xmin)+error;
			xmax=min(p0[2][j][i].xmax,p0[2][j-1][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j-1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j-1][i1].zmax)-error;
			if(xmin<xmax&&zmin<zmax)
			{
			    norm=1;
				r1=min(dminy(&p0[2][j-1][i1],&p0[2][j][i]),p0[2][j][i].ymin-p0[2][j-1][i1].ymin);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(j>=2)
	{
		for(i1=initialy[0];i1<=finaly[0];i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j-2][i1].xmin)+error;
			xmax=min(p0[2][j][i].xmax,p0[2][j-2][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j-2][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j-2][i1].zmax)-error;
			if(xmin<xmax&&zmin<zmax)
			{
			    norm=1;
				r1=dminy(&p0[2][j-2][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(i>=1)
	{
		for(i1=initialx;i1<=finalx;i1++)
		{
			xmin=max(p0[2][j][i].xmin,p0[2][j][i1].xmin)+error;
			xmax=min(p0[2][j][i].xmax,p0[2][j][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j][i1].zmax)-error;
			if(xmin<xmax&&zmin<zmax)
			{
				norm=1;
				r1=dminy(&p0[2][j][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}
    if(norm==0)
		r=0;
	return(r);
};


double EllipseGen::dminy(struct parameter0 *p1,struct parameter0 *p2)
{
	struct coordinate point[21][21];
	double xmin,xmax,zmin,zmax;
	double r=0,r1,r2;
	int i,j;
	int norm1=0,norm2=0;

	xmin=max((*p1).xmin,(*p2).xmin)+error;
	xmax=min((*p1).xmax,(*p2).xmax)-error;
	zmin=max((*p1).zmin,(*p2).zmin)+error;
	zmax=min((*p1).zmax,(*p2).zmax)-error;
	r1=(*p2).ymin-(*p1).ymin;
	 r2=r1;
    
	for(j=0;j<21;j++)
		for(i=0;i<21;i++)
		{
			point[j][i].x=xmin+i*(xmax-xmin)/20;
			point[j][i].z=zmin+j*(zmax-zmin)/20;
		}

	for(j=0;j<21;j++)
	{
		for(i=0;i<21;i++)
		{
			if(judgey(p1,p2,point[j][i].x,point[j][i].z))
			{
	    		r=disy(p1,p2,point[j][i].x,point[j][i].z);
					
				if(r<r1)
				{
		    	   r1=r;
				}
			}
			else
				continue;
					
		}
	}


	return(r1);

};

double EllipseGen::disy(struct parameter0 *p1,struct parameter0 *p2,double x,double z)
{
	double r;
	double A1,B1,C1,D1,E1,F1;
	double A2,B2,C2,D2,E2,F2;
	A1=(pow((*p1).coe[3],2)-(*p1).coe[1]*(*p1).coe[0])/pow((*p1).coe[1],2);
	A2=(pow((*p2).coe[3],2)-(*p2).coe[1]*(*p2).coe[0])/pow((*p2).coe[1],2);
	B1=(pow((*p1).coe[4],2)-(*p1).coe[2]*(*p1).coe[1])/pow((*p1).coe[1],2);
    B2=(pow((*p2).coe[4],2)-(*p2).coe[2]*(*p2).coe[1])/pow((*p2).coe[1],2);
	C1=((*p1).coe[3]*(*p1).coe[4]-(*p1).coe[1]*(*p1).coe[5])/pow((*p1).coe[1],2);
	C2=((*p2).coe[3]*(*p2).coe[4]-(*p2).coe[1]*(*p2).coe[5])/pow((*p2).coe[1],2);
	D1=(*p1).coe[6]/(*p1).coe[1];
	D2=(*p2).coe[6]/(*p2).coe[1];
	E1=-(*p1).coe[3]/(*p1).coe[1];
	E2=-(*p2).coe[3]/(*p2).coe[1];
    F1=-(*p1).coe[4]/(*p1).coe[1];
	F2=-(*p2).coe[4]/(*p2).coe[1];

	r=-sqrt(A2*pow(x-(*p2).center_x,2)+B2*pow(z-(*p2).center_z,2)+2*C2*(x-(*p2).center_x)*(z-(*p2).center_z)+D2)
		        +E2*(x-(*p2).center_x)+F2*(z-(*p2).center_z)+(*p2).center_y
	   -sqrt(A1*pow(x-(*p1).center_x,2)+B1*pow(z-(*p1).center_z,2)+2*C1*(x-(*p1).center_x)*(z-(*p1).center_z)+D1)
	            -E1*(x-(*p1).center_x)-F1*(z-(*p1).center_z)-(*p1).center_y;
	return(r);
}

int EllipseGen::judgey(struct parameter0 *p1,struct parameter0 *p2,double x,double z)
{
	double A1,B1,C1,D1;
	double A2,B2,C2,D2;
	int norm;

	A1=(*p1).coe[1]*(*p1).coe[0]-pow((*p1).coe[3],2);
	B1=(*p1).coe[1]*(*p1).coe[2]-pow((*p1).coe[4],2);
	C1=(*p1).coe[1]*(*p1).coe[5]-(*p1).coe[3]*(*p1).coe[4];
	D1=(*p1).coe[6]*(*p1).coe[1];

	A2=(*p2).coe[1]*(*p2).coe[0]-pow((*p2).coe[3],2);
	B2=(*p2).coe[1]*(*p2).coe[2]-pow((*p2).coe[4],2);
	C2=(*p2).coe[1]*(*p2).coe[5]-(*p2).coe[3]*(*p2).coe[4];
	D2=(*p2).coe[6]*(*p2).coe[1];

	if(A1*pow(x-(*p1).center_x,2)+B1*pow(z-(*p1).center_z,2)+C1*2*(x-(*p1).center_x)*(z-(*p1).center_z)<=D1
	    &&A2*pow(x-(*p2).center_x,2)+B2*pow(z-(*p2).center_z,2)+C2*2*(x-(*p2).center_x)*(z-(*p2).center_z)<=D2)
		norm=1;
	else
		norm=0;

	return(norm);
}

/*******************************************************************************************************/
/*                             y定位完毕                             */
 /******************************************************************************************************/


  /**********************************************************************/
  /*                    x的第二次定位                                */
 /**********************************************************************/

double EllipseGen::secondsettlex(struct parameter0 ***p0,int k,int j,int i)
{
	double r1,r;
	int i1,j1;
	double ymin,ymax,zmin,zmax;
    int norm=0;

	r=min(p0[2][j][i].xmin-p0[2][j][i-1].xmin,p0[2][j][i].xmax-p0[2][j][i-1].xmax);

	for(j1=initialwz[1];j1<=finalwz[1];j1++)
		for(i1=initialkz[1];i1<=finalkz[1];i1++)
		{
			ymin=max(p0[2][j][i].ymin,p0[1][j1][i1].ymin)+error;
			ymax=min(p0[2][j][i].ymax,p0[1][j1][i1].ymax)-error;
			zmin=max(p0[2][j][i].zmin,p0[1][j1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[1][j1][i1].zmax)-error;
			if(ymin<ymax&&zmin<zmax)
			{
				norm=1;
				r1=dminx(&p0[1][j1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}

	if(k>1)
	{
		for(j1=initialwz[0];j1<=finalwz[0];j1++)
		for(i1=initialkz[0];i1<=finalkz[0];i1++)
		{
			ymin=max(p0[2][j][i].ymin,p0[0][j1][i1].ymin)+error;
			ymax=min(p0[2][j][i].ymax,p0[0][j1][i1].ymax)-error;
			zmin=max(p0[2][j][i].zmin,p0[0][j1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[0][j1][i1].zmax)-error;
			if(ymin<ymax&&zmin<zmax)
			{
			    norm=1;
				r1=dminx(&p0[0][j1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}


	if(j>=1)
	{
		for(i1=initialy[1];i1<=finaly[1];i1++)
		{
			ymin=max(p0[2][j][i].ymin,p0[2][j-1][i1].ymin)+error;
			ymax=min(p0[2][j][i].ymax,p0[2][j-1][i1].ymax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j-1][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j-1][i1].zmax)-error;
			if(ymin<ymax&&zmin<zmax)
			{
			    norm=1;
				r1=dminx(&p0[2][j-1][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(j>=2)
	{
		for(i1=initialy[0];i1<=finaly[0];i1++)
		{
			ymin=max(p0[2][j][i].xmin,p0[2][j-2][i1].xmin)+error;
			ymax=min(p0[2][j][i].xmax,p0[2][j-2][i1].xmax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j-2][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j-2][i1].zmax)-error;
			if(ymin<ymax&&zmin<zmax)
			{
			    norm=1;
				r1=dminx(&p0[2][j-2][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}

	if(i>=1)
	{
		for(i1=initialx;i1<=finalx;i1++)
		{
			ymin=max(p0[2][j][i].ymin,p0[2][j][i1].ymin)+error;
			ymax=min(p0[2][j][i].ymax,p0[2][j][i1].ymax)-error;
			zmin=max(p0[2][j][i].zmin,p0[2][j][i1].zmin)+error;
			zmax=min(p0[2][j][i].zmax,p0[2][j][i1].zmax)-error;
			if(ymin<ymax&&zmin<zmax)
			{
			    norm=1;
				r1=dminx(&p0[2][j][i1],&p0[2][j][i]);
				if(r1<r)
					r=r1;
			}
		}
	}
 //   printf("r %lf\n",r);
	if(norm==0)
		r=0;
	return(r);
};


double EllipseGen::dminx(struct parameter0 *p1,struct parameter0 *p2)
{
	struct coordinate point[21][21];
	double ymin,ymax,zmin,zmax;
	double r=0,r1,r2;
	int i,j;
	int norm1=0,norm2=0;

	ymin=max((*p1).ymin,(*p2).ymin)+error;
	ymax=min((*p1).ymax,(*p2).ymax)-error;
	zmin=max((*p1).zmin,(*p2).zmin)+error;
	zmax=min((*p1).zmax,(*p2).zmax)-error;
	r1=(*p2).xmin-(*p1).xmin;
	 r2=r1;
    
	for(j=0;j<21;j++)
		for(i=0;i<21;i++)
		{
			point[j][i].y=ymin+i*(ymax-ymin)/20;
			point[j][i].z=zmin+j*(zmax-zmin)/20;
		}

	for(j=0;j<21;j++)
	{
		for(i=0;i<21;i++)
		{
			if(judgex(p1,p2,point[j][i].y,point[j][i].z))
			{
	    		r=disx(p1,p2,point[j][i].y,point[j][i].z);
					
				if(r<r1)
				{
		    	   r1=r;
				}
			}
			else
				continue;
					
		}
	}


	return(r1);

};

double EllipseGen::disx(struct parameter0 *p1,struct parameter0 *p2,double y,double z)
{
	double r;
	double A1,B1,C1,D1,E1,F1;
	double A2,B2,C2,D2,E2,F2;
	A1=(pow((*p1).coe[3],2)-(*p1).coe[1]*(*p1).coe[0])/pow((*p1).coe[0],2);
	A2=(pow((*p2).coe[3],2)-(*p2).coe[1]*(*p2).coe[0])/pow((*p2).coe[0],2);
	B1=(pow((*p1).coe[5],2)-(*p1).coe[2]*(*p1).coe[0])/pow((*p1).coe[0],2);
    B2=(pow((*p2).coe[5],2)-(*p2).coe[2]*(*p2).coe[0])/pow((*p2).coe[0],2);
	C1=((*p1).coe[3]*(*p1).coe[5]-(*p1).coe[0]*(*p1).coe[4])/pow((*p1).coe[0],2);
	C2=((*p2).coe[3]*(*p2).coe[5]-(*p2).coe[0]*(*p2).coe[4])/pow((*p2).coe[0],2);
	D1=(*p1).coe[6]/(*p1).coe[0];
	D2=(*p2).coe[6]/(*p2).coe[0];
	E1=-(*p1).coe[3]/(*p1).coe[0];
	E2=-(*p2).coe[3]/(*p2).coe[0];
    F1=-(*p1).coe[5]/(*p1).coe[0];
	F2=-(*p2).coe[5]/(*p2).coe[0];

	r=-sqrt(A2*pow(y-(*p2).center_y,2)+B2*pow(z-(*p2).center_z,2)+2*C2*(y-(*p2).center_y)*(z-(*p2).center_z)+D2)
		        +E2*(y-(*p2).center_y)+F2*(z-(*p2).center_z)+(*p2).center_x
	   -sqrt(A1*pow(y-(*p1).center_y,2)+B1*pow(z-(*p1).center_z,2)+2*C1*(y-(*p1).center_y)*(z-(*p1).center_z)+D1)
	            -E1*(y-(*p1).center_y)-F1*(z-(*p1).center_z)-(*p1).center_x;
	return(r);
}

int EllipseGen::judgex(struct parameter0 *p1,struct parameter0 *p2,double y,double z)
{
	double A1,B1,C1,D1;
	double A2,B2,C2,D2;
	int norm;

	A1=(*p1).coe[1]*(*p1).coe[0]-pow((*p1).coe[3],2);
	B1=(*p1).coe[0]*(*p1).coe[2]-pow((*p1).coe[5],2);
	C1=(*p1).coe[0]*(*p1).coe[4]-(*p1).coe[3]*(*p1).coe[5];
	D1=(*p1).coe[6]*(*p1).coe[0];

	A2=(*p2).coe[1]*(*p2).coe[0]-pow((*p2).coe[3],2);
	B2=(*p2).coe[0]*(*p2).coe[2]-pow((*p2).coe[5],2);
	C2=(*p2).coe[0]*(*p2).coe[4]-(*p2).coe[3]*(*p2).coe[5];
	D2=(*p2).coe[6]*(*p2).coe[0];

	if(A1*pow(y-(*p1).center_y,2)+B1*pow(z-(*p1).center_z,2)+C1*2*(y-(*p1).center_y)*(z-(*p1).center_z)<=D1
	    &&A2*pow(y-(*p2).center_y,2)+B2*pow(z-(*p2).center_z,2)+C2*2*(y-(*p2).center_y)*(z-(*p2).center_z)<=D2)
		norm=1;
	else
		norm=0;

	return(norm);
}
/*******************************************************************************************************/
/*                             x定位完毕                             */
 /******************************************************************************************************/



int EllipseGen::inserty(struct parameter0 *p0,struct parameter0 *p)
{
	struct parameter0 p1,p2;
	int i1,k1=N2;

	if((*p0).xmin<p[0].xmin)
		k1=-1;

    for(i1=0;i1<N2;i1++)
	{
		if((*p0).xmin>=p[i1].xmin&&(p[i1+1].a==0||i1==N2-1))
		{
            k1=i1;
			break;
		}
		else if((*p0).xmin>=p[i1].xmin&&(*p0).xmin<p[i1+1].xmin)
		{
			k1=i1;
			break;
		}
	}
	
	if(k1>=N2-1)
		return 0;

	replace(&p1,p0);
    for(i1=k1+1;i1<N2;i1++)
	{
		replace(&p2,&p[i1]);
		replace(&p[i1],&p1);
		replace(&p1,&p2);
	}

	return 1;
    
};




int EllipseGen::insertz(struct parameter0 *p0,struct parameter0 **p)
{
	int i1,j1,k1,k2;
	struct parameter0 p1,p2;

	if((*p0).ymin<total_ymin[1][0])
		k1=0;

	for(j1=0;j1<N1;j1++)
	{
		if((*p0).ymin>=total_ymin[1][j1]&&(total_ymax[1][j1+1]==0||j1==N1-1))
		{
            k1=j1;
			break;
		}
		else if((*p0).ymin>=total_ymin[1][j1]&&(*p0).ymin<total_ymin[1][j1+1])
		{
			k1=j1;
			break;
		}
	}

	
	if((*p0).xmin<p[k1][0].xmin)
		k2=-1;

    for(i1=0;i1<N2;i1++)
	{
		if((*p0).xmin>=p[k1][i1].xmin&&(p[k1][i1+1].al==0||i1==N2-1))
		{
            k2=i1;//printf("k1 k2 p[k1][k2] %d %d %lf %lf\n",k1,k2,p[k1][k2].xmin,p[k1][k2+1].xmin);
			break;
			
		}
		else if((*p0).xmin>=p[k1][i1].xmin&&(*p0).xmin<p[k1][i1+1].xmin)
		{
			k2=i1;//printf("k1 k2 p[k1][k2] %d %d %lf %lf\n",k1,k2,p[k1][k2].xmin,p[k1][k2+1].xmin);
			break;
		}
	}
    

	if(k1>=N1-1||k2>=N2-1)
		return 0;

	replace(&p1,p0);
    for(i1=k2+1;i1<N2;i1++)
	{
		replace(&p2,&p[k1][i1]);
		replace(&p[k1][i1],&p1);
		replace(&p1,&p2);
	}

   // printf("k1 k2 %d %d\n",k1,k2);

	if((*p0).ymin<total_ymin[1][k1])
		total_ymin[1][k1]=(*p0).ymin;
	if((*p0).ymax>total_ymax[1][k1])
		total_ymax[1][k1]=(*p0).ymax;                      



	for(i1=k2+1;i1<N2;i1++)
	{
	
		for(j1=0;j1<N1;j1++)
		{
			if(p[j1][i1].al!=0)
			{
				total_xmin[1][i1]=p[j1][i1].xmin;
				total_xmax[1][i1]=p[j1][i1].xmax;
			    break;
			}
		
		}
		for(;j1<N1;j1++)
		{
			if(p[j1][i1].xmin<total_xmin[1][i1]&&p[j1][i1].al!=0)
				total_xmin[1][i1]=p[j1][i1].xmin;
			if(p[j1][i1].xmax>total_xmax[1][i1]&&p[j1][i1].al!=0)
				total_xmax[1][i1]=p[j1][i1].xmax;
		}
	}

	return 1;
};




void EllipseGen::replace(struct parameter0 *p1,struct parameter0 *p2)
{
	int k1;
	
		  {
			  (*p1).A=(*p2).A;
			  (*p1).a=(*p2).a;
			  (*p1).al=(*p2).al;
			  (*p1).alpha1=(*p2).alpha1;
			  (*p1).alpha2=(*p2).alpha2;
			  (*p1).alpha3=(*p2).alpha3;
			  (*p1).B=(*p2).B;
			  (*p1).b=(*p2).b;
			  (*p1).bata1=(*p2).bata1;
			  (*p1).bata2=(*p2).bata2;
			  (*p1).bata3=(*p2).bata3;
			  (*p1).bl=(*p2).bl;
			  (*p1).C=(*p2).C;
			  (*p1).c=(*p2).c;
			  (*p1).center_x=(*p2).center_x;
			  (*p1).center_y=(*p2).center_y;
			  (*p1).center_z=(*p2).center_z;
			  (*p1).cl=(*p2).cl;

			  for(k1=0;k1<=6;k1++)
			  {
	    		  (*p1).coe[k1]=(*p2).coe[k1];
			  }

			  (*p1).gama1=(*p2).gama1;
			  (*p1).gama2=(*p2).gama2;
			  (*p1).gama3=(*p2).gama3;
			  (*p1).xmax=(*p2).xmax;
			  (*p1).ymax=(*p2).ymax;
			  (*p1).xmin=(*p2).xmin;
			  (*p1).ymin=(*p2).ymin;
			  (*p1).zmin=(*p2).zmin;
			  (*p1).zmax=(*p2).zmax;

		  }

};





void EllipseGen::parametergeneration1(struct parameter0 *p1,int f_e_c_division,int distribution_sign1)
{
	double X1,X2,X3;
	double miu1=0,miu2=0,gama1=0,gama2=0;
	double amin,amax,bmin,cmin;
	double error1=pow(10.0,-7);
	double *iaseed=NULL,*ibseed=NULL,*icseed=NULL;

	if(f_e_c_division==0)
	{
		amin=famin;
		amax=famax;
		iaseed=&f_aseed;
	}
	else if(f_e_c_division==1)
	{
		amin=eamin;
		amax=eamax;
		bmin=ebmin;
		cmin=ecmin;
		iaseed=&e_aseed;
		ibseed=&e_bseed;
		icseed=&e_cseed;
	}
	else
	{
		amin=camin;
		amax=camax;
		bmin=cbmin;
		cmin=cc;
		iaseed=&c_aseed;
		ibseed=&c_bseed;
	}
	(*p1).a=unifrnd(amin,amax,iaseed);
    (*p1).al=(*p1).a+h*(*p1).a;

	/*if(f_e_c_division==1)
       (*p1).b=unifrnd(bmin,(*p1).a,ibseed);
	else if(f_e_c_division==2)
        (*p1).b=bmin;
	else if(f_e_c_division==0)
		(*p1).b=fb;

    (*p1).bl=(*p1).b+h*(*p1).b;
	if(f_e_c_division==0)
	{
		(*p1).c=(*p1).b;
	}
	else if(f_e_c_division==1)
	{
		(*p1).c=unifrnd(cmin,(*p1).b,icseed);
	}
	else
	{
		(*p1).c=cc;
	}
	(*p1).cl=(*p1).c+h*(*p1).c;*/

	(*p1).b=bmin;
	(*p1).c=cmin;
	(*p1).bl=(*p1).b+h*(*p1).b;
	(*p1).cl=(*p1).c+h*(*p1).c;

	if(distribution_sign1==0)
	{
        gama1=unifrnd(gamamin,gamamax,&ig1seed);
		miu1=unifrnd(miumin,miumax,&im1seed);
	}
	else
	{
		
		double t1=gasdev(&ig1seed);
		double t2=gasdev(&im1seed);
		gama1=sqrt(sitagama)*t1+miugama;
		miu1=sqrt(sitamiu)*t2+miumiu;

		//将gama1变换做0-pi内的值，miu1变换做0-2pi内的值
		if(gama1<0)
		{
			gama1=2*pi-fmod(-gama1,2*pi);
		}
		if(miu1<0)
		{
			miu1=2*pi-fmod(-miu1,2*pi);
		}
		gama1=fmod(gama1,2*pi);
		if(gama1>pi&&gama1<2*pi)
		{
			gama1=gama1-pi;
		}
		miu1=fmod(miu1,2*pi);
	}

	(*p1).gama1=cos(gama1);
	(*p1).alpha1=cos(miu1)*sin(gama1);
    (*p1).bata1=sin(miu1)*sin(gama1);


	if(gama1>=0&&gama1<=pi/2)
	{
		if(f_e_c_division==1||f_e_c_division==2)
	    	gama2=unifrnd(pi/2-gama1,pi/2+gama1,&ig2seed);
		else
			gama2=pi/2;
        (*p1).gama2=cos(gama2);
	}
	else if(gama1>pi/2&&gama1<=pi)
	{
		if(f_e_c_division==1||f_e_c_division==2)
    		gama2=unifrnd(gama1-pi/2,3*pi/2-gama1,&ig2seed);
		else
			gama2=pi/2;
		(*p1).gama2=cos(gama2);
	}
   

	if(fabs(fabs((*p1).gama1)-1)<error)
	{
		if(f_e_c_division==1||f_e_c_division==2)
	    	miu2=unifrnd(0,2*pi,&im2seed);
		else
			miu2=pi;
		(*p1).alpha2=cos(miu2);
        (*p1).bata2=sin(miu2);   
	}
	else if(fabs(fabs((*p1).gama2)-1)<error)
	{
		(*p1).alpha2=0;
		(*p1).bata2=0; 

	}
    else
	{
		double cosmiu2;
		X1=(1-pow((*p1).gama1,2))*sqrt(1-pow((*p1).gama2,2));
		X2=-(*p1).gama1*(*p1).gama2*(*p1).alpha1;
		X3=pow((*p1).alpha1,2)*pow((*p1).gama1,2)*pow((*p1).gama2,2)
			+(1-pow((*p1).alpha1,2))*(1-pow((*p1).gama1,2))*(1-pow((*p1).gama2,2))-(1-pow((*p1).gama1,2))*pow((*p1).gama1,2);
		cosmiu2=(X2+sqrt(fabs(X3)))/X1;

        (*p1).alpha2=cosmiu2*sin(gama2);

		if(fabs((*p1).bata1)<error)
			(*p1).bata2=sqrt(fabs(1-pow((*p1).alpha2,2)-pow((*p1).gama2,2)));
		else
            (*p1).bata2=-((*p1).gama1*(*p1).gama2+(*p1).alpha1*(*p1).alpha2)/(*p1).bata1;	
	}            

    //printf("miu1 gama1 gama2 %lf %lf %lf\n",miu1,gama1,gama2);
    
	(*p1).gama3=(*p1).alpha1*(*p1).bata2-(*p1).alpha2*(*p1).bata1;
    (*p1).bata3=-(*p1).gama2*(*p1).alpha1+(*p1).alpha2*(*p1).gama1;
  	(*p1).alpha3=-(*p1).gama1*(*p1).bata2+(*p1).bata1*(*p1).gama2;

    
   	
    (*p1).A=pow((*p1).al,2)*pow((*p1).alpha1,2)+pow((*p1).bl,2)*pow((*p1).alpha2,2)
		                                                        +pow((*p1).cl,2)*pow((*p1).alpha3,2);
	(*p1).B=pow((*p1).al,2)*pow((*p1).bata1,2)+pow((*p1).bl,2)*pow((*p1).bata2,2)
		                                                        +pow((*p1).cl,2)*pow((*p1).bata3,2);
	(*p1).C=pow((*p1).al,2)*pow((*p1).gama1,2)+pow((*p1).bl,2)*pow((*p1).gama2,2)
		                                                        +pow((*p1).cl,2)*pow((*p1).gama3,2);
};



void EllipseGen::parametergeneration2(struct parameter0 *p0)
{
	double A1,A2,A3,B1,B2,B3,C1,C2,C3;
	double a,b,c;
	double F;
	A1=((*p0).bata3)*((*p0).gama2)-((*p0).bata2)*((*p0).gama3);
    A2=((*p0).bata1)*((*p0).gama3)-((*p0).bata3)*((*p0).gama1);
	A3=((*p0).bata2)*((*p0).gama1)-((*p0).bata1)*((*p0).gama2);
	B1=((*p0).alpha2)*((*p0).gama3)-((*p0).alpha3)*((*p0).gama2);
    B2=((*p0).alpha3)*((*p0).gama1)-((*p0).alpha1)*((*p0).gama3);
	B3=((*p0).alpha1)*((*p0).gama2)-((*p0).alpha2)*((*p0).gama1);
	C1=((*p0).alpha3)*((*p0).bata2)-((*p0).alpha2)*((*p0).bata3);
	C2=((*p0).alpha1)*((*p0).bata3)-((*p0).alpha3)*((*p0).bata1);
	C3=((*p0).alpha2)*((*p0).bata1)-((*p0).alpha1)*((*p0).bata2);
	F=B1*((*p0).bata1)+B2*((*p0).bata2)+B3*((*p0).bata3);

    a=(*p0).al;
	b=(*p0).bl;
	c=(*p0).cl;


	(*p0).coe[0]=pow(a,2)*pow(b,2)*pow(A3,2)+pow(a,2)*pow(c,2)*pow(A2,2)+pow(b,2)*pow(c,2)*pow(A1,2);
	(*p0).coe[1]=pow(a,2)*pow(b,2)*pow(B3,2)+pow(a,2)*pow(c,2)*pow(B2,2)+pow(b,2)*pow(c,2)*pow(B1,2);
	(*p0).coe[2]=pow(a,2)*pow(b,2)*pow(C3,2)+pow(a,2)*pow(c,2)*pow(C2,2)+pow(b,2)*pow(c,2)*pow(C1,2);
	(*p0).coe[3]=pow(a,2)*pow(b,2)*A3*B3+pow(a,2)*pow(c,2)*A2*B2+pow(b,2)*pow(c,2)*A1*B1;
	(*p0).coe[4]=pow(a,2)*pow(b,2)*C3*B3+pow(a,2)*pow(c,2)*C2*B2+pow(b,2)*pow(c,2)*C1*B1;
	(*p0).coe[5]=pow(a,2)*pow(b,2)*A3*C3+pow(a,2)*pow(c,2)*A2*C2+pow(b,2)*pow(c,2)*A1*C1;
	(*p0).coe[6]=pow(a,2)*pow(b,2)*pow(c,2)*pow(F,2);
};


double EllipseGen::min(double x0,double x1)
{
  return (x0<x1)?x0:x1;;
};



double EllipseGen::max(double x0,double x1)
{
  return (x0>x1)?x0:x1;;
};



double EllipseGen::unifrnd(double c,double d,double* iseed)
{
	double r;
	r=recursion(iseed)/M;
	r=c+(d-c)*r;
	return r;
};


double EllipseGen::recursion(double* iseed)
{
  double random; 
  random=fmod(lamda*(*iseed),M);
  *iseed=random;
  return random;
}
;

double EllipseGen::gasdev(double *iseed)
{
	static int iset;
	static double gset;
	double t,v1,v2,fac,r;
	if(iset==0)
	{
		do
		{
			v1=2.0*unifrnd(0.0,1.0,iseed)-1.0;
			v2=2.0*unifrnd(0.0,1.0,iseed)-1.0;
			r=v1*v1+v2*v2;
		}while((r>=1.0)||(r==0.0));
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		t=v2*fac;
		iset=1;
		return t;
	
	}
	else
	{
		t=gset;
		iset=0;
		return t;
	}
}

double EllipseGen::normfunction(double t0,double t1,double x)
{
  double r;
  r=exp(-pow(x-t0,2)/(2*pow(t1,2)))/(t1*sqrt(2*pi));
  return r;
}


double EllipseGen::expfunction(double t,double x)
{
	double r;
	r=t*exp(-t*x);
	return r;
}

char* EllipseGen::get_line(ifstream& input_stream, char* rline){
   //     char rline[250];     //读入的一行
        input_stream.getline(rline,250);
        //跳过注释行
        while( !input_stream.eof() && rline[0] == '%' )
                input_stream.getline(rline,250);
  //      yout << "rline: " << rline << endl;
        return rline;
};


double EllipseGen::volume(double x1,double y1,double z1,struct parameter *p0)
{
	double xmin,ymin,zmin,xmax,ymax,zmax;
	double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
	double W[4],W1;
	int norm,i;
	double a,b,c;
	double A,B,C;
	a=(*p0).a;
	b=(*p0).b;
	c=(*p0).c;
	A=(*p0).A;
	B=(*p0).B;
	C=(*p0).C;
	xmin=x1-sqrt(A);
    xmax=x1+sqrt(A);
	ymin=y1-sqrt(B);
	ymax=y1+sqrt(B);
	zmin=z1-sqrt(C);
	zmax=z1+sqrt(C);
	Xmin=xmin-x1;
	Xmax=xmax-x1;
	Ymin=ymin-y1;
	Ymax=ymax-y1;
	Zmin=zmin-z1;
	Zmax=zmax-z1;

	for(i=0;i<=3;i++)
		W[i]=0;
	

	if(xmin<leftbar&&xmax>leftbar)
	{
		double X=leftbar-x1;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		norm=0;
	}
	else if(xmax>rightbar&&xmin<rightbar)
	{
		double X=rightbar-x1;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		W[0]=4*pi*a*b*c/3-W[0];
		norm=1;
	}
	 if(ymin<backbar&&ymax>backbar)
	{
		double Y=backbar-y1;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		norm=2;
	}
	else if(ymax>formbar&&ymin<formbar)
	{
		double Y=formbar-y1;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		W[1]=4*pi*a*b*c/3-W[1];
		norm=3;
	}
	 if(zmin<downbar&&zmax>downbar)
	{
		double Z=downbar-z1;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		norm=4;
	}
	else if(zmax>upbar&&zmin<upbar)
	{
		double Z=upbar-z1;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		W[2]=4*pi*a*b*c/3-W[2];
		norm=5;
	}
	if(xmin>leftbar&&xmax<rightbar&&ymin>backbar&&ymax<formbar&&zmin>downbar&&zmax<upbar)
	{
		W[3]=4*pi*a*b*c/3;
		norm=6;
	}

	W1=(rightbar-leftbar)*(formbar-backbar)*(upbar-downbar);
	if(norm!=6)
	{
		for(i=0;i<=2;i++)
		{
			if(W[i]!=0&&W[i]<W1)
				W1=W[i];
		}
	}
	else
		W1=W[3];
//	printf("W1 %lf\n",W1);

	return(W1);
	
}


double EllipseGen::volume(double x1,double y1,double z1,struct basicparameter *p0)
{
	double xmin,ymin,zmin,xmax,ymax,zmax;
	double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
	double W[4],W1;
	int norm,i;
	double a,b,c;
	double A,B,C;
	a=(*p0).a;
	b=(*p0).b;
	c=(*p0).c;
	A=(*p0).A;
	B=(*p0).B;
	C=(*p0).C;
	xmin=x1-sqrt(A);
    xmax=x1+sqrt(A);
	ymin=y1-sqrt(B);
	ymax=y1+sqrt(B);
	zmin=z1-sqrt(C);
	zmax=z1+sqrt(C);
	Xmin=xmin-x1;
	Xmax=xmax-x1;
	Ymin=ymin-y1;
	Ymax=ymax-y1;
	Zmin=zmin-z1;
	Zmax=zmax-z1;

	for(i=0;i<=3;i++)
		W[i]=0;
	

	if(xmin<leftbar&&xmax>leftbar)
	{
		double X=leftbar-x1;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		norm=0;
	}
	else if(xmax>rightbar&&xmin<rightbar)
	{
		double X=rightbar-x1;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		W[0]=4*pi*a*b*c/3-W[0];
		norm=1;
	}
	 if(ymin<backbar&&ymax>backbar)
	{
		double Y=backbar-y1;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		norm=2;
	}
	else if(ymax>formbar&&ymin<formbar)
	{
		double Y=formbar-y1;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		W[1]=4*pi*a*b*c/3-W[1];
		norm=3;
	}
	 if(zmin<downbar&&zmax>downbar)
	{
		double Z=downbar-z1;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		norm=4;
	}
	else if(zmax>upbar&&zmin<upbar)
	{
		double Z=upbar-z1;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		W[2]=4*pi*a*b*c/3-W[2];
		norm=5;
	}
	if(xmin>leftbar&&xmax<rightbar&&ymin>backbar&&ymax<formbar&&zmin>downbar&&zmax<upbar)
	{
		W[3]=4*pi*a*b*c/3;
		norm=6;
	}

	W1=(rightbar-leftbar)*(formbar-backbar)*(upbar-downbar);
	if(norm!=6)
	{
		for(i=0;i<=2;i++)
		{
			if(W[i]!=0&&W[i]<W1)
				W1=W[i];
		}
	}
	else
		W1=W[3];
//	printf("W1 %lf\n",W1);

	return(W1);
	
}

//生成随机数的初值
double EllipseGen::seed(void)
{
	//用于随机生成的起始时间
	return ((double)(rand())/RAND_MAX)*M;
}

