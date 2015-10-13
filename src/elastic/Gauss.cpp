//--------------------------------------------
//高斯类实现；
//Gauss.cpp
//---------------------------------------------

#include"Gauss.h"

int Gauss::Generate_gauss(int type){
	if(type==2){
               vector<Node> tem_gauss;
		       vector<double> tem_wight;
			   if(precision==2){
				   tem_gauss.push_back(Node(0.25,0.25,0.25));
				   tem_wight.push_back(1.0);
			   }
			   else if(precision==3){
				   double a=0.58541020;
				   double b=0.13819660;
				   tem_gauss.push_back(Node(a,b,b));
				   tem_gauss.push_back(Node(b,a,b));
				   tem_gauss.push_back(Node(b,b,a));
				   tem_gauss.push_back(Node(b,b,b));

				   tem_wight.push_back(0.25);
				    tem_wight.push_back(0.25);
				    tem_wight.push_back(0.25);
				    tem_wight.push_back(0.25);
			   }
			   gauss=tem_gauss;
               wight=tem_wight;
	}
 else if(type==3){
     if(precision==2){
		     vector<double> tem_gaussz;
		     vector<double> tem_wightz;
			 //-------------------------------------
	 double a=0.7745966692;
			tem_gaussz.push_back(-a);
			tem_gaussz.push_back(0);
			tem_gaussz.push_back(a);
			//--------------------------------------
			a=0.5555555556;
	 double b=0.8888888889;
			tem_wightz.push_back(a);
			tem_wightz.push_back(b);
			tem_wightz.push_back(a);
			//--------------------------------------
			 vector<Vector2D> tem_gaussxy;
			 vector<double> tem_wightxy;
			 tem_gaussxy.push_back(Vector2D(0.5,0.5));
			 tem_gaussxy.push_back(Vector2D(0.5,0));
			 tem_gaussxy.push_back(Vector2D(0,0.5));
            //---------------------------------------
			 tem_wightxy.push_back(0.3333333333);
			 tem_wightxy.push_back(0.3333333333);
			 tem_wightxy.push_back(0.3333333333);
            //---------------------------------------
		for(int i=0;i<int(tem_gaussz.size());i++)
	       for(int j=0;j<int(tem_gaussxy.size());j++){
				wight.push_back(tem_wightz[i]*tem_wightxy[j]);
				Node n;
				n.x=tem_gaussxy[j].x;
				n.y=tem_gaussxy[j].y;
				n.z=tem_gaussz[i];
		     gauss.push_back(n);
		}
	}
	}

return 1;
}
//=====================================================================
