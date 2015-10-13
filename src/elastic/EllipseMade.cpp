//===========================================================================
// EllipseMade.cpp
// �����������������Ա�������ų�����
// Member Functions in a class of ellipse mading with random
//===========================================================================
#include "EllipseMade.h"

const double	 pi=3.1415926535897;
const double	 epsi=0.01;
const double	 epsilon=0.01;
const int N_times=100000;


//---------------------------------------------------------------------------
void EllipseMade::ellip_generation(ifstream &infile,string output_file,string data_file,double mod)
{
	struct  point_type
		{
			double	x;
			double	y;
			double	z;
		};

	double	x,y,z;
	double	x1,y1,z1;
	double	A,B,C;
	double	a_max,a_min,b_min,c_min;		// the ellipse the max a and the min a and the min b,a is half of the long axis,b is the half of the short axis.		                
	double  sita,phi;
	double  alpha1,beta1,alpha2;
	double	sign;
	double  d,f;								// the distance between two ellipses center
	double  vol_ratio;
	int		times;
	int		k1,K,l1,L;

	istringstream istr0(Get_Line(infile));		//�����������������
        istr0	>> space.x0 >> space.x >> space.y0 >> space.y >> space.z0 >> space.z;	

	istringstream istr2(Get_Line(infile));		//�����ж�����Ϣ	
		istr2 >> a_min >> a_max >> b_min >> c_min;
	if (a_min>a_max)
	{
		hout << "�������������Сֵ�������ֵ��" << endl;
		exit(0);
	}

	istringstream istr6(Get_Line(infile));		//����ٷֱ�
	istr6 >> vol_ratio;
//@
	//	int in;	
	//	ifstream o("num.dat");
	//	o>>in;
	//	o.close();
	//if(in>10) in=in-10;
	//else if(in>5) in=in-5;
	//vol_ratio=vol_ratio*in;

	//-------------------------------------------------------------
	//�������ɼ�����
	double delt_h = 2*sqrt(pow(c_min,2)*(pow(a_max,2)-pow(a_max-(a_min+c_min)*epsilon,2))/pow(a_max,2))/pi; 
	double cell_volume = (space.x-space.x0)*(space.y-space.y0)*(space.z-space.z0);
	double ellip_volume = 0.0;
	ellip_ratio = 0.0;
	times=0;

	//����������ɵ���ʼʱ��
	srand((unsigned int)time(0));

	do
	{
		//-------------------------------------------------------------
		//��������
		struct elliparam ell_temp;

		ell_temp.x=((double)(rand())/RAND_MAX)*(space.x-space.x0)+space.x0;
		ell_temp.y=((double)(rand())/RAND_MAX)*(space.y-space.y0)+space.y0;
		ell_temp.z=((double)(rand())/RAND_MAX)*(space.z-space.z0)+space.z0;

		ell_temp.a=((double)(rand())/RAND_MAX)*(a_max-a_min)+a_min;
		if(!(b_min==0&&c_min==0))
		{
			ell_temp.b=((double)(rand())/RAND_MAX)*(ell_temp.a-b_min)+b_min;
			ell_temp.c=((double)(rand())/RAND_MAX)*(ell_temp.a-c_min)+c_min;
		}
		else
		{
			ell_temp.b=ell_temp.a;
			ell_temp.c=ell_temp.a;
		}

		alpha1 =((double)(rand())/RAND_MAX)*pi;
		if(alpha1>pi/2.0)
		{
			beta1  =((double)(rand())/RAND_MAX)*(2*(pi-alpha1))+(alpha1-pi/2.0);
		}
		else
		{
			beta1  =((double)(rand())/RAND_MAX)*(2*alpha1)+(pi/2.0-alpha1);
		}
		ell_temp.alpha1 =cos(alpha1);																// r1 from 0 to pi 
		ell_temp.beta1  =cos(beta1);																	// r2 from (pi/2-r1) to (pi/2+r1)
		ell_temp.gamma1 =(int)pow(-1.0,(int)fmod(rand(),2.0)+1)*sqrt(1-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	// choose one in both, minus value or positive value

		if(alpha1>pi/2.0)
		{
			alpha2  =((double)(rand())/RAND_MAX)*(2*(pi-alpha1))+(alpha1-pi/2.0);
		}
		else
		{
			alpha2  =((double)(rand())/RAND_MAX)*(2*alpha1)+(pi/2.0-alpha1);
		}
		ell_temp.alpha2=cos(alpha2);

		A=1+pow(ell_temp.beta1/ell_temp.gamma1,2);
		B=2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
		C=pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;

		ell_temp.beta2	 =(-B+(int)pow(-1.0,(int)fmod(rand(),2.0)+1)*sqrt(pow(B,2)-4*A*C))/(2.0*A);
		ell_temp.gamma2 =-(ell_temp.beta1/ell_temp.gamma1)*ell_temp.beta2-(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1);
		
		sign=(ell_temp.alpha1*ell_temp.beta2)/fabs(ell_temp.alpha1*ell_temp.beta2);
		ell_temp.alpha3 =sign*sqrt(1-pow(ell_temp.alpha1,2)-pow(ell_temp.alpha2,2));
		ell_temp.beta3	 =-(ell_temp.alpha1*ell_temp.beta1+ell_temp.alpha2*ell_temp.beta2)/ell_temp.alpha3;
		ell_temp.gamma3 =-(ell_temp.alpha1*ell_temp.gamma1+ell_temp.alpha2*ell_temp.gamma2)/ell_temp.alpha3;

		ell_temp.a=(1+epsilon)*ell_temp.a;
		ell_temp.b=(1+epsilon)*ell_temp.b;
		ell_temp.c=(1+epsilon)*ell_temp.c;
		//---------------------------------------------------------------------------
		//check this ellipse
		//---------------------------------------------------------------------------
		//whether this ellipse intersect with other ellipses
		k1=(int)(sqrt(pow(ell_temp.a,2)+pow(ell_temp.b,2))/delt_h);
		K=4*(k1+1);
		sita=pi/(K/2.0);

		//(���һ�γ���080414)���ǲ���߽��ཻ�����
		//for(int j=1;j<=K/2;j++)
		//{
		//	l1=(int)(sqrt(pow(ell_temp.a*sin(j*sita),2)+pow(ell_temp.b*sin(j*sita),2))/delt_h);
		//	L=4*(l1+1);
		//	phi=2*pi/L;

		//	for(int m=1;m<=L;m++)
		//	{
		//		x=ell_temp.a*sin(j*sita)*cos(m*phi);
		//		y=ell_temp.b*sin(j*sita)*sin(m*phi);
		//		z=ell_temp.c*cos(j*sita);

		//		x1=ell_temp.x+x*ell_temp.alpha1+y*ell_temp.alpha2+z*ell_temp.alpha3;
		//		y1=ell_temp.y+x*ell_temp.beta1+y*ell_temp.beta2+z*ell_temp.beta3;
		//		z1=ell_temp.z+x*ell_temp.gamma1+y*ell_temp.gamma2+z*ell_temp.gamma3;

		//		x=x1;
		//		y=y1;
		//		z=z1;

		//		if(x1<=space.x0+0.04||x1>=space.x-0.04||y1<=space.y0+0.04||y1>=space.y-0.04||z1<=space.z0+0.04||z1>=space.z-0.04)
		//		{
		//			times=times+1;
		//			goto gen_again;
		//		}
		//	}
		//}

		for(int i=0;i<int(ellip.size());i++)
		{
			//probably estimate
			d=sqrt(pow(ell_temp.x-ellip[i].x,2)+pow(ell_temp.y-ellip[i].y,2)+pow(ell_temp.z-ellip[i].z,2));

			if(d>=ell_temp.a+ellip[i].a+2*epsi)
			{
				goto gene;
			}
			else if((d<ell_temp.b+ellip[i].b+2*epsi*(ellip[i].b/ellip[i].a))&&(d<ell_temp.c+ellip[i].c+2*epsi*(ellip[i].c/ellip[i].a)))
			{
				times=times+1;
				goto gen_again;
			}
			else 
			{
				//exactly compute
				for(int j=1;j<=K/2;j++)
				{
					l1=(int)(sqrt(pow(ell_temp.a*sin(j*sita),2)+pow(ell_temp.b*sin(j*sita),2))/delt_h);
					L=4*(l1+1);
					phi=2*pi/L;
		
					for(int m=1;m<=L;m++)
					{
						x=ell_temp.a*sin(j*sita)*cos(m*phi);
						y=ell_temp.b*sin(j*sita)*sin(m*phi);
						z=ell_temp.c*cos(j*sita);

						x1=ell_temp.x+x*ell_temp.alpha1+y*ell_temp.alpha2+z*ell_temp.alpha3;
						y1=ell_temp.y+x*ell_temp.beta1+y*ell_temp.beta2+z*ell_temp.beta3;
						z1=ell_temp.z+x*ell_temp.gamma1+y*ell_temp.gamma2+z*ell_temp.gamma3;

						x=x1;
						y=y1;
						z=z1;

						x=x-ellip[i].x;
						y=y-ellip[i].y;
						z=z-ellip[i].z;
			
						x1=x*ellip[i].alpha1+y*ellip[i].beta1+z*ellip[i].gamma1;
						y1=x*ellip[i].alpha2+y*ellip[i].beta2+z*ellip[i].gamma2;
						z1=x*ellip[i].alpha3+y*ellip[i].beta3+z*ellip[i].gamma3;
			
						f=pow(x1,2)/pow(ellip[i].a,2)+pow(y1,2)/pow(ellip[i].b,2)+pow(z1,2)/pow(ellip[i].c,2)-1.0;

						if(f<0.0)
						{
							times=times+1;
							goto gen_again;
						}
					}
				}
			}
gene: ;
		}
		//---------------------------------------------------------------------
		//�������㲢������������
		times=0;
		ellip.push_back(ell_temp);
		//---------------------------------------------------------------------
		//���������
		ellip_volume = ellip_volume + ellipse_volume(ell_temp,epsilon);
		ellip_ratio = ellip_volume/cell_volume;
gen_again:	 ;
	}while(times<=N_times&&ellip_ratio<vol_ratio);

	//---------------------------------------------------------------------
	//��������
	for(int i=0;i<int(ellip.size());i++)
	{
		ellip[i].a=ellip[i].a/(1+epsilon);
		ellip[i].b=ellip[i].b/(1+epsilon);
		ellip[i].c=ellip[i].c/(1+epsilon);
	}

	//---------------------------------------------------------------------
	//��ͼ
//	if(mod!=0) draw_tecplot("generate_ellipse.dat",delt_h);

	//---------------------------------------------------------------------
	//�������
	if(mod!=0) output_EllipseData(output_file);
	
	//---------------------------------------------------------------------
	//���뵽���������ļ�
	output_Datafile(data_file);

}
//---------------------------------------------------------------------------
//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
string EllipseMade::Get_Line(ifstream &input_stream)const
{  
	string str;
	//�ȶ���һ����Ϣ
	getline(input_stream,str);
	//����ע����
	while(!input_stream.eof() && str.substr(0,1)=="%")
	{
		getline(input_stream,str);
	}
	return str;
}
//---------------------------------------------------------------------------
//���������ڵ����е����
double EllipseMade::ellipse_volume(const struct elliparam &ellip,const double &epsilon)
{
	double xmin,ymin,zmin,xmax,ymax,zmax;
	double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
	double W[4],W1;
	int norm;
	double a,b,c;
	double A,B,C;

	a=ellip.a/(1+epsilon);
	b=ellip.b/(1+epsilon);
	c=ellip.c/(1+epsilon);

    A=pow(a,2)*pow(ellip.alpha1,2)+pow(b,2)*pow(ellip.alpha2,2)+pow(c,2)*pow(ellip.alpha3,2);
	B=pow(a,2)*pow(ellip.beta1,2)+pow(b,2)*pow(ellip.beta2,2)+pow(c,2)*pow(ellip.beta3,2);
	C=pow(a,2)*pow(ellip.gamma1,2)+pow(b,2)*pow(ellip.gamma2,2)+pow(c,2)*pow(ellip.gamma3,2);

	xmin=ellip.x-sqrt(A);
    xmax=ellip.x+sqrt(A);
	ymin=ellip.y-sqrt(B);
	ymax=ellip.y+sqrt(B);
	zmin=ellip.z-sqrt(C);
	zmax=ellip.z+sqrt(C);
	Xmin=xmin-ellip.x;
	Xmax=xmax-ellip.x;
	Ymin=ymin-ellip.y;
	Ymax=ymax-ellip.y;
	Zmin=zmin-ellip.z;
	Zmax=zmax-ellip.z;

	for(int i=0;i<=3;i++)	W[i]=0;
	
	if(xmin<space.x0&&xmax>space.x0)
	{
		double X=space.x0-ellip.x;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		norm=0;
	}
	else if(xmax>space.x&&xmin<space.x)
	{
		double X=space.x-ellip.x;
		W[0]=pi*a*b*c*(Xmax-X-(pow(Xmax,3)-pow(X,3))/(3*A))/sqrt(A);
		W[0]=4*pi*a*b*c/3-W[0];
		norm=1;
	}
	 if(ymin<space.y0&&ymax>space.y0)
	{
		double Y=space.y0-ellip.y;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		norm=2;
	}
	else if(ymax>space.y&&ymin<space.y)
	{
		double Y=space.y-ellip.y;
		W[1]=pi*a*b*c*(Ymax-Y-(pow(Ymax,3)-pow(Y,3))/(3*B))/sqrt(B);
		W[1]=4*pi*a*b*c/3-W[1];
		norm=3;
	}
	 if(zmin<space.z0&&zmax>space.z0)
	{
		double Z=space.z0-ellip.z;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		norm=4;
	}
	else if(zmax>space.z&&zmin<space.z)
	{
		double Z=space.z-ellip.z;
		W[2]=pi*a*b*c*(Zmax-Z-(pow(Zmax,3)-pow(Z,3))/(3*C))/sqrt(C);
		W[2]=4*pi*a*b*c/3-W[2];
		norm=5;
	}
	if(xmin>space.x0&&xmax<space.x&&ymin>space.y0&&ymax<space.y&&zmin>space.z0&&zmax<space.z)
	{
		W[3]=4*pi*a*b*c/3;
		norm=6;
	}

	W1=(space.x-space.x0)*(space.y-space.y0)*(space.z-space.z0);
	if(norm!=6)
	{
		for(int i=0;i<=2;i++)
		{
			if(W[i]!=0&&W[i]<W1)	W1=W[i];
		}
	}
	else
	{
		W1=W[3];
	}

	return(W1);
}

//---------------------------------------------------------------------------
//��ͼ
void EllipseMade::draw_tecplot(const string &outfile, const double &delt_h)const
{
	ofstream out(outfile.c_str());
	for(int i=0; i<int(ellip.size()); i++)
	{
		out << "title = generate_ellipse" << endl;
		out << "variables = x, y, z" << endl;
		out << "zone t= \" line \" " << endl;

		double x,y,z;
		double x1,y1,z1;
		double k1=(int)(sqrt(pow(ellip[i].a,2)+pow(ellip[i].b,2))/delt_h);
		double K=4*(k1+1);
		double sita=pi/(K/2.0);
		for(int j=1;j<=K/2;j++)
		{
			double l1=(int)(sqrt(pow(ellip[i].a*sin(j*sita),2)+pow(ellip[i].b*sin(j*sita),2))/delt_h);
			double L=4*(l1+1);
			double phi=2*pi/L;
			
			for(int m=1;m<=L;m++)
			{
				x=ellip[i].a*sin(j*sita)*cos(m*phi);
				y=ellip[i].b*sin(j*sita)*sin(m*phi);
				z=ellip[i].c*cos(j*sita);

				x1=ellip[i].x+x*ellip[i].alpha1+y*ellip[i].alpha2+z*ellip[i].alpha3;
				y1=ellip[i].y+x*ellip[i].beta1+y*ellip[i].beta2+z*ellip[i].beta3;
				z1=ellip[i].z+x*ellip[i].gamma1+y*ellip[i].gamma2+z*ellip[i].gamma3;

				x=x1;
				y=y1;
				z=z1;

				out << x << "  " << y << "  " << z << endl;
			}

			x=ellip[i].a*sin(j*sita)*cos(phi);
			y=ellip[i].b*sin(j*sita)*sin(phi);
			z=ellip[i].c*cos(j*sita);

			x1=ellip[i].x+x*ellip[i].alpha1+y*ellip[i].alpha2+z*ellip[i].alpha3;
			y1=ellip[i].y+x*ellip[i].beta1+y*ellip[i].beta2+z*ellip[i].beta3;
			z1=ellip[i].z+x*ellip[i].gamma1+y*ellip[i].gamma2+z*ellip[i].gamma3;

			x=x1;
			y=y1;
			z=z1;

			out << x << "  " << y << "  " << z << endl;
		}
	}
	out << "title = generate_ellipse" << endl;
	out << "variables = x, y, z" << endl;
	out << "zone n=8, e=1, f=fepoint, et=brick" << endl;
	out << space.x0 <<"  " << space.y0 <<"  " << space.z0 <<"  " << endl;
	out << space.x <<"  " << space.y0 <<"  " << space.z0 <<"  " << endl;
	out << space.x <<"  " << space.y <<"  " << space.z0 <<"  " << endl;
	out << space.x0 <<"  " << space.y <<"  " << space.z0 <<"  " << endl;
	out << space.x0 <<"  " << space.y0 <<"  " << space.z<<"  "  << endl;
	out << space.x <<"  " << space.y0 <<"  " << space.z <<"  " << endl;
	out << space.x <<"  " << space.y <<"  " << space.z <<"  " << endl;
	out << space.x0 <<"  " << space.y <<"  " << space.z <<"  " << endl;
	out << "1 2 3 4 5 6 7 8" << endl;
	out.close();

}
//---------------------------------------------------------------------
//�������
void EllipseMade::output_EllipseData(const string &outfile)const
{
	ofstream out(outfile.c_str());
	out <<"%�����ڿ����������������" << endl;
	out << int(ellip.size()) << " " << ellip_ratio <<endl;
	out <<"%������ʼ������" <<endl;
	out <<space.x0 <<" " <<space.y0 <<" "<<space.z0 <<endl;
	out <<"%���������" <<endl;
	out <<space.x <<" " <<space.y <<" "<<space.z <<endl;
	out <<"%�������15������:�����������ꡢ����ߴ硢��λ��" << endl;
	for(int i=0; i<int(ellip.size ()); i++)
	{
		out	<< ellip[i].x << " " << ellip[i].y << " " << ellip[i].z << " "
				<< ellip[i].a << " " << ellip[i].b << " " << ellip[i].c << " "
				<< ellip[i].alpha1 << " " << ellip[i].alpha2 << " " << ellip[i].alpha3 << " "
				<< ellip[i].beta1 << " " << ellip[i].beta2 << " " << ellip[i].beta3 << " "
				<< ellip[i].gamma1 << " " << ellip[i].gamma2 << " " << ellip[i].gamma3 << " "
				<< endl;
	}
	out.close();
}
//---------------------------------------------------------------------
//���뵽���������ļ�
void EllipseMade::output_Datafile(const string data_file)const
{
	ofstream out(data_file.c_str(),ios::app);
	out <<"%�����ڿ����������������" << endl;
	out << int(ellip.size()) << " " << ellip_ratio <<endl;
	out.close();
}

//===========================================================================
