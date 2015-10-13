//===========================================================================
// SolveEqu.cpp
// 解线性方程组函数
//===========================================================================
#include "SolveEqu.h"

//---------------------------------------------------------------------------
//紧缩存储刚度矩阵（动态申请空间+izig函数）
int SolveEqu::Gen_izig(const vector<Node> &nodes_vec, const vector<Element> &elements_vec, int* &Iz, int* &Ig, const int &LR)
{
	//---------------------------------------------------------------------------
	//动态申请存储空间
	int N		=	(int)nodes_vec.size();	//节点个数
	int NTEL = (int)elements_vec.size();	//单元个数

	int EPNum = 0;		//单元中节点的总个数
	for (int i=0; i < NTEL; i++)
	{
		EPNum = EPNum+int(elements_vec[i].nodes_id.size());
	}

	//Iem: 单元材料号
	int *Iem = new int[NTEL];
	assert(Iem);

	//Ien: 单元节点号
	int *Ien = new int[EPNum];
	assert(Ien);

	//转换单元信息
	int epn=0;
	for(int i=0; i<NTEL; i++ )
	{
		Iem[i] = int(elements_vec[i].nodes_id.size());			//单元节点个数
		for (int j=0; j<Iem[i]; j++)
		{
			Ien[epn+j] = elements_vec[i].nodes_id[j]+1;
		}
		epn=epn+Iem[i];
	}
	if (epn!=EPNum)
	{
		hout << "给Ien赋单元节点编号时，总个数不对！" << endl;
		return 0;
	}

	//紧缩存储刚度矩阵
	//以下程序izig主要输出数组iz和ig，是变带宽存储
	//数组iz的长度N，ig的长度是IB1，IB1是程序izig内部一个数据，小于LR
	izig(N, NTEL, Ien, Iem, Iz, Ig, LR);

	//释放存储空间
	delete Iem;
	delete Ien;

	return 1;
}
//---------------------------------------------------------------------------
//紧缩存储刚度矩阵
void SolveEqu::izig(	const int &N, const int &NEL, int* &IEN, int* &IEM, int *&IZ, int *&IG,const int &LR	)
{
	int IC[80];

	for(int I=0; I<N; I++)
	{
		IZ[I]=0;	//IZ[K-2]中记录与第K个点相关的数据在IG中存储的开始位置
	}
	int IB1=0;		//累加器
	int IB2;			

	for(int K=2; K <= N; K++)
	{
		famnss(	K,IEN,IEM,&IB2,IC,NEL);
		if( IB2 != 0 )
		{
			int L2 = IZ[K-2]+IB2; //IB2记录与K相关的节点（同在一个单元中）中，整体编号比K小的节点的个数（无重复）
			if( L2 >= LR )	//LR是一个大数，在本程序中取MAXIG*N
			{
				hout << " warning: L2 >= LR, in IZIG.cpp " << endl;
				cout << " warning: L2 >= LR, in IZIG.cpp " << endl;
				exit(0);
			}
			else
			{
				for(int L=1; L <= IB2; L++)
				{
					IG[IB1+L-1] = IC[L]; //添加到IG，IC中存放与K号节点相关的节点（同在一个单元中）中，整体编号比K小的节点的编号
				}
				IB1=IB1+IB2;	//计数器累加
				IZ[K-1] = L2;  //L2==IB1,IZ[K-2]和IZ[K-1]中分别记录与第K个点相关的数据在IG中存储的开始位置和结束位置
			}
		}
		else
		{
			IZ[K-1] = IZ[K-2]; //没有相关且编号小的节点，开始和结束是一点
		}
	}
}
//---------------------------------------------------------------------------
//在函数izig中被调用
void SolveEqu::famnss(int INOD, int *IEN, int *IEM, int *IA, int *KN, int NEL)
{
	int NOD1,NOD2;

	int K1=0;
	*IA=0;

	for( int J=0; J < NEL; J++ )	//遍历单元，NEL是单元个数
	{
		for( int J1=1; J1 <= IEM[J]; J1++ )	//遍历单元中的节点，IEM[J]是第J个单元的节点个数
		{
			if( IEN[K1+J1-1] == 0 ) //IEN[K1+J1-1]是第K1+J1-1个节点的整体编号；K1是累加数；因为IEN数组从零开始，所以J1要减1
			{									  //整体编号从1开始，所以等于0是错的
				hout << "出错！！,出现把四面体单元当作三棱柱单元处理节点的错误！" << endl;
				exit(0);
			}
			if( IEN[K1+J1-1] == INOD )   //INOD 是调用这个函数时那个节点的编号，等于K
			{												//找到K号节点所在单元
				for( int J2=1; J2 <= IEM[J]; J2++ )
				{
					NOD2=IEN[K1+J2-1];
					if(NOD2 < INOD && NOD2 != 0) //在这个单元中的其他节点的整体编号小于K时
					{
						if(*IA != 0)		//在KN中查找是否有编号重复的节点，*IA记录与K相关的节点（同在一个单元中）中，整体编号比K小的节点的个数（无重复）
						{						//KN中放的是与K号节点相关的节点（同在一个单元中）中，整体编号比K小的节点的编号
							for( int J3=1; J3 <= *IA; J3++ )
							{
								if(NOD2 == KN[J3])
								{
									goto label1;
								}
							}
						}
						*IA=*IA+1;
						KN[*IA]=NOD2;  //插入KN数组
					}
label1: ;
				}
				goto label2;
			}
		}
label2:		K1=K1+IEM[J];  //节点个数累加
	}


	if( *IA >= 2 )	//KN从小到大排序，KN中放的是与K号节点相关的节点（同在一个单元中）中，整体编号比K小的
	{
		int NN;
		for( int L=2; L <= *IA; L++ )
		{
			NN=L-1;
			NOD2=KN[L];
			for( int JN=1; JN <= NN; JN++ )
			{
				NOD1=KN[JN];
				if( NOD2 < NOD1 )
				{
					for( int JJN=NN; JJN >= JN; JJN-- )
					{
						KN[JJN+1]=KN[JJN];
					}
					KN[JN]=NOD2;
					break;
				}
			}
		}
	}
}

//---------------------------------------------------------------------------
//处理位移边界条件
void SolveEqu::tacp(	int M14, int N0, int *IZ, int *IG, int *IP, double *VP, double *F, double *AK )
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * *     
//     CALL TOCP2.																							*     
//     CALLED BY MAIN.																					*     
//     PURPOSE:																								*     
//     TREAT THREE_DIMENSION CONSTRAINT POINTS WITH ZERO.	*     
//     M14_THE NUMBER OF CONSTRAINT INFORMATIONS.					*     
//     N0_LAST NODE OF CURRENT STRATIFICATION.								*     
//     IZ_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.				*     
//     IG_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.				*     
//     IP_CONSTRAINT POINT AND TYPE.													*     
//     VP_CONTRAINT VALUES.																	*     
//     F_THE EXTERNAL FORCES.																	*     
//     AK_STIFFNESS MATRIX.																		*     
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
// LKt   u    v    w
//   7     0    0    0   及其它这种 固定边值条件
//	 6    空 0    0    空表示这一维是自由边界  其它0处，可以是0或非零值 但必须有值
//   5     0   空 0    
//	 4   空 空 0
//   3     0    0  空
//   2   空 0   空
//   1    0   空 空

{
	int IA,IH;

	//初始化对象
	Mt=0;	LKt=0;
	Jt1=Jt2=Jt3=Jt4=Jt5=Jt6=Jt7=Jt8=Jt9=0;

	if( M14 == 0 )	goto L90;
	for( int K=1; K <= M14; K++ )
	{
		IA=2*(K-1);                                                      
		IH=3*(K-1);                                                       
		Mt=IP[IA];
		if( Mt > N0 )	goto L80;
		LKt=IP[IA+1];
		switch( LKt )
		{
		case 1: goto L10;
		case 2: goto L20;
		case 3: goto L30;
		case 4: goto L40;
		case 5: goto L50;
		case 6: goto L60;
		case 7: goto L70;
		case 8: goto L80;
		default: break;
		}
L10:	if( VP[IH] != 0.0 )	goto L80;
L11:	Jt1=1;                                     
		Jt2=3;                                                              
		Jt5=1;                                                              
		Jt6=7;                                                              
		Jt9=3; 
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L20:	if( VP[IH+1] != 0.0 )	goto L80;
L21:	Jt1=4;                                                                                         
		Jt2=6;                                                              
		Jt5=2;                                                             
		Jt6=8;                                                              
		Jt9=3;   
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L30:	if( VP[IH] == 0.0 && VP[IH+1] == 0.0 )	goto L32;
		if( VP[IH] != 0.0 )	goto L31;                   
		LKt=1;                                                              
		goto L11;
L31:	if( VP[IH+1] != 0.0 )	goto L80;                   
		LKt=2;                                                              
		goto L21;                                                          
L32:	Jt1=1;                                                              
		Jt2=3;                                                              
		Jt3=4;                                                            
		Jt4=6;                                                              
		Jt5=1;                                                              
		Jt6=7;                                                              
		Jt7=2;                                                              
		Jt8=8;                                                              
		Jt9=3;                                                              
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L40:	if( VP[IH+2] != 0.0 )	goto L80;                   
L41:	Jt1=7;                                                              
		Jt2=9;                                                              
		Jt5=3;                                                              
		Jt6=9;                                                              
		Jt9=3;                                                              
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L50:  if( VP[IH] == 0.0 && VP[IH+2] == 0.0 )	goto L52;
		if( VP[IH] != 0.0 )	goto L51;                                         
		LKt=1;                                                            
		goto L11;                                                          
L51:  if( VP[IH+2] != 0.0 )	goto L80;                            
		LKt=4;                                                              
		goto L41;                                                          
L52:	Jt1=1;                                                             
		Jt2=3;                                                              
		Jt3=7;                                                              
		Jt4=9;                                                              
		Jt5=1;                                                              
		Jt6=7;                                                              
		Jt7=3;                                                              
		Jt8=9;                                                              
		Jt9=3;                                                              
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L60:  if( VP[IH+1] == 0.0 && VP[IH+2] == 0.0 )	goto L62;
		if( VP[IH+1] != 0.0 )	goto L61;                                         
		LKt=2;                                                            
		goto L21;                                                          
L61:  if( VP[IH+2] != 0.0 )	goto L80;                            
		LKt=4;                                                              
		goto L41;                                                                                                                
L62:	Jt1=4;                                                              
		Jt2=6;                                                              
		Jt3=7;                                                              
		Jt4=9;                                                              
		Jt5=2;                                                              
		Jt6=8;                                                              
		Jt7=3;                                                              
		Jt8=9;                                                              
		Jt9=3;                                                              
		tocp2(N0,IZ,IG,F,AK);
		goto L80;
L70:	if( VP[IH] == 0.0 && VP[IH+1] == 0.0 && VP[IH+2] == 0.0 )	goto L32;
		if( VP[IH] == 0.0 && VP[IH+1] == 0.0 )	goto L76;
		if( VP[IH] == 0.0 && VP[IH+2] == 0.0 )	goto L75;
		if( VP[IH+1] == 0.0 && VP[IH+2] == 0.0 )	goto L74;
		if( VP[IH] == 0.0 )	goto L73; 
		if( VP[IH+1] == 0.0 )	goto L72;
		if( VP[IH+2] == 0.0 )	goto L71; 
		goto L80;   
L71:	LKt=4;                                                             
		goto L41;
L72:	LKt=2;                                                             
		goto L21;
L73:	LKt=1;                                                             
		goto L11;
L74:	LKt=6;                                                             
		goto L62;
L75:	LKt=5;                                                             
		goto L52;
L76:	LKt=3;                                                             
		goto L32;
L80:	;
	}
L90:	;
}
//---------------------------------------------------------------------------
//在函数tacp中被调用
void SolveEqu::tocp2( int N0, int *IZ, int *IG, double *F, double *AK )
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     
//     CALLED BY TACP.																		*     
//     PURPOSE:																					*     
//     TREAT ONE CONSTRINT POINT WITH ZERO.						*     
//     N0_LAST NODE OF CURRENT STRATIFICATION.					*     
//     IZ_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.	*     
//     IG_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.	*     
//     F_THE EXTERNAL FORCES.														*     
//     AK_STIFFNESS MATRIX.															*     
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     
{
		int IC,M1;

		switch( LKt )
		{
		case 1: goto L10;
		case 2: goto L20;
		case 3: goto L30;
		case 4: goto L40;
		case 5: goto L50;
		case 6: goto L60;
		case 7: goto L70;
		default: break;
		}                                
L10:	F[3*(Mt-1)]=0.0;                                            
		goto L80;                                                          
L20:	F[3*(Mt-1)+1]=0.0;                                           
		goto L80;                                                          
L30: F[3*(Mt-1)]=0.0;                                         
		F[3*(Mt-1)+1]=0.0;                                        
		goto L80;                                                          
L40:	F[3*(Mt-1)+2]=0.0;                                     
		goto L80;                                                        
L50:	F[3*(Mt-1)]=0.0;                                       
		F[3*(Mt-1)+2]=0.0;                                      
		goto L80;                                                        
L60:	F[3*(Mt-1)+1]=0.0;                                      
		F[3*(Mt-1)+2]=0.0;                                     
		goto L80;                                                
L70:	F[3*(Mt-1)]=0.0;
		F[3*(Mt-1)+1]=0.0;
		F[3*(Mt-1)+2]=0.0;
		goto L80;      
L80:	if (Mt == 1)
		{
			IC = 0;
		}
		else
		{
			IC=IZ[Mt-1]-IZ[Mt-2];
		}
		if( IC == 0 ) goto L110;
		int KK;
		for( int L=1; L <= IC; L++ )
		{
			for( int J=Jt1; J <= Jt2; J++ )
			{
				KK=6*(Mt-1)+9*(IZ[Mt-2]+L-1)+J;
				AK[KK-1]=0.0;
			}
			if( LKt == 3 || LKt == 5 || LKt == 6 || LKt == 7 ) goto L90;
			goto L105;
L90:		for( int J=Jt3; J <= Jt4; J++ )
			{
				KK=6*(Mt-1)+9*(IZ[Mt-2]+L-1)+J;
				AK[KK-1]=0.0;
			}
			if( LKt != 7 )		goto L105;
			for( int J=7; J <= 9; J++ )
			{
				KK=6*(Mt-1)+9*(IZ[Mt-2]+L-1)+J;
				AK[KK-1]=0.0;
			}
L105:	;
		}

L110:M1=Mt+1;
		for( int I=M1; I <= N0; I++ )
		{
			IC=IZ[I-1]-IZ[I-2];
			if( IC == 0 )	goto L140;
			for( int L=1; L <= IC; L++ )
			{
				if( IG[IZ[I-2]+L-1] != Mt )	goto L135;
				for( int J=Jt5; J <= Jt6; J=J+Jt9 )
				{
					KK=6*(I-1)+9*(IZ[I-2]+L-1)+J;                                      
					AK[KK-1]=0.0;
				}
				if( LKt == 3 || LKt == 5 || LKt == 6 || LKt == 7 ) goto L120;
				goto L135;
L120:		for( int J=Jt7; J <= Jt8; J=J+Jt9 )
				{
					KK=6*(I-1)+9*(IZ[I-2]+L-1)+J;                                      
					AK[KK-1]=0.0;
				}
				if( LKt != 7 )	goto L135; 
				for( int J=3; J <= 9; J=J+3 )
				{
					KK=6*(I-1)+9*(IZ[I-2]+L-1)+J;                                      
					AK[KK-1]=0.0;
				}
L135: ;	}
L140: ;
		}
		switch( LKt )
		{
		case 1: goto L150;
		case 2: goto L151;
		case 3: goto L152;
		case 4: goto L153;
		case 5: goto L154;
		case 6: goto L155;
		case 7: goto L156;
		default: break;
		}
L150:KK=9*IZ[Mt-1]+6*(Mt-1)+1;                                              
		AK[KK-1]=1.0;
		AK[KK]=0.0;
		AK[KK+2]=0.0;
		goto L160;                                                         
L151:KK=9*IZ[Mt-1]+6*(Mt-1)+2;                                              
		AK[KK-1]=0.0;
		AK[KK]=1.0;
		AK[KK+2]=0.0;
		goto L160;                                                        
L152:KK=9*IZ[Mt-1]+6*(Mt-1)+1;                                              
		AK[KK-1]=1.0;
		AK[KK]=0.0;
		AK[KK+1]=1.0;
		AK[KK+2]=0.0;
		AK[KK+3]=0.0;
		goto L160;                                                       
L153:KK=9*IZ[Mt-1]+6*(Mt-1)+4;                                              
		AK[KK-1]=0.0;
		AK[KK]=0.0;
		AK[KK+1]=1.0;
		goto L160;                                                        
L154:KK=9*IZ[Mt-1]+6*(Mt-1)+1;                                              
		AK[KK-1]=1.0;
		AK[KK]=0.0;
		AK[KK+2]=0.0;
		AK[KK+3]=0.0;
		AK[KK+4]=1.0;
		goto L160;                                                        
L155:KK=9*IZ[Mt-1]+6*(Mt-1)+2;                                              
		AK[KK-1]=0.0;
		AK[KK]=1.0;
		AK[KK+1]=0.0;
		AK[KK+2]=0.0;
		AK[KK+3]=1.0;
		goto L160;                                                       
L156:KK=9*IZ[Mt-1]+6*(Mt-1)+1;                                              
		AK[KK-1]=1.0;
		AK[KK]=0.0;
		AK[KK+1]=1.0;
		AK[KK+2]=0.0;
		AK[KK+3]=0.0;
		AK[KK+4]=1.0;	
L160:	;
}

//---------------------------------------------------------------------------
// 求解线性方程组函数
void SolveEqu::sol( double EX, int M14, int N, int *IZ, int *IG, int *IP, double *VP, double *A, 
				   double *B, double *P, double *R, double *S, double *V, double *X )
				   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     
				   //     CALL TCNOZ,MABVM.																					*     
				   //     CALLED BY MAIN.																							*     
				   //     PURPOSE:																										*     
				   //     THIS ROUTINE IS THE SOLUTION OF LINEAR EQUATIONS BY THE		*     
				   //     CONJUGATE GRADIENT METHOD.																*     
				   //     EX_CONTRAL ERROR.																					*     
				   //     KM_ITERATIVE KM TIMES PRINT MIDDLE RESULTS.								*     																				*     
				   //     M14_THE NUMBER OF CONSTRAINT POINTS.											*     
				   //     N__LAST NODE OF CURRENT STRATIFICATION.										*     
				   //     N1_2*N OR 3*N.																								*     
				   //     IZ_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.						*     
				   //     IG_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.						*     
				   //     IP_CONSTRAINT POINT AND TYPE.															*     
				   //     VP-CONSTRAINT VALUES.																			*     
				   //     A--STIFFNESS MATRIX.																				*     
				   //     B,P,R,S,V_WORKING STORAGE.																	*     
				   //     X_RESULT STORAGE.																					*     
				   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
{
	int N1,NPR,NTI,NUP,K,KK,ISF;
	double Rk,R0,RR0,APP,AK,RRk,BK;
	N1=3*N;
	NPR=100;                  
	R0=N1;   
	Rk=0.0;                                               
	NTI=240*N1;                               
	NUP=240*N1; 
	K=0; 
	//**********************************************************************
	//                                                                     
	//---- PROCESS CONSTRAINTS FOR MATRIX A,VECTOR X,B                    
	//**********************************************************************
	if( M14 != 0 ) tcnoz( 1,1,M14,N,IP,VP,X,V );
L10:	;
	KK=0;
	//----	CALCULATE PRODUCT A*X0=> AP  
	mabvm( N,N1,IZ,IG,A,X,V );
	//---- CALCULATE R0=P0=B-A*X0
	for( int I=0; I < N1; I++ )
	{
		R[I]=B[I]-V[I];
		P[I]=R[I];
	}
	//----	PROCESS CONSTRAINTS FOR P0,R0 ABOUT SPECIAL TEMPERATURES
	if( M14 != 0 ) tcnoz(0,1,M14,N,IP,VP,P,R);
	//---------------------------------------------------------------------
	//----	COMPUTE R(K)*R(K)
	//---------------------------------------------------------------------
	RR0=0.0;
	for( int I=0; I < N1; I++ )
	{
		RR0=RR0+R[I]*R[I];
	}
	if( K == 0 ) R0=sqrt(RR0);
	//--------------------------------------------------------------------- 
	//---- ENTER TO ITERATE                                                
	//---------------------------------------------------------------------
L50:	K=K+1;
	KK=KK+1;
	//
	//---- CALCULATE AP=A*P=>S
	//
	mabvm( N,N1,IZ,IG,A,P,S );
	//
	//----	PROCESS CONSTRAINTS FOR AP ABOUT SPECIAL TEMPERATURES       
	//                                                               
	//      IF(M14.NE.0) CALL TCNOZ(0,0,M14,N,IP,VP,S,P)  
	//                                                              
	//---- CALCULATE INNER PRODUCT (AP,P)                          
	//
	APP=0.0;
	for( int I=0; I < N1; I++ )
	{
		APP=APP+S[I]*P[I];
	}                                    
	AK=-RR0/APP;                                     
	//
	//---- CALCULATE X(K)=X(K-1)-AK*P(K)                  
	//---- CALCULATE R(K)=R(K-1)+AK*AP(K)                
	//
	RRk=0.0;
	for( int I=0; I < N1; I++ )
	{
		X[I]=X[I]-AK*P[I];
		R[I]=R[I]+AK*S[I];
		RRk=RRk+R[I]*R[I];
	}
	Rk=sqrt(RRk);
	if( fabs(Rk) >= EX*fabs(R0) )
	{
		if( K >= NUP )
		{
			//------	ITERATION DIVERGE
			ISF=1;
			goto L240;
		}
		else
		{
			BK=RRk/RR0;
			RR0=RRk;
			for( int I=0; I < N1; I++ )
			{
				P[I]=R[I]+BK*P[I];
			}
			if( KK == NTI ) goto L10;
			if( KK < NTI ) goto L50;
		}
	}
	else
	{
		if( (int)fmod(K*1.0,NPR*1.0) != 0 )	ISF=3;
	}
L240:	;
}

//---------------------------------------------------------------------------
//在函数sol中被调用
void SolveEqu::tcnoz( int KG, int KG1, int M14, int N, int *IP, double *VP, double *W, double *G )
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *    
//     CALLED BY SOL.																			*     
//     PURPOSE:																						*     
//     TREAT CONSTRAINT POINTS WITH NON_ZERO.						*     
//     KG_MARK.																						*     
//     ZERO SENDS W, OR VP SENDS W.	    										*     
//     KG1_WHEN KG1 EQUALS ZERO,ZERO SENDS G.                      * 
//     M14_THE NUMBER OF CONSTRAINT POINTS.							*     
//     N__LAST NODE OF CURRENT STRATIFICATION.						*     
//     IP_CONSTRAINT POINT AND TYPE.											*     
//     VP_CONSTRAINT VALUES.															*     
//     W_RESULT STORAGE.																	*     
//     G_RESULT STORAGE.																	*     
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * *     
// LKt   u    v    w
//   7     0    0    0   及其它这种 固定边值条件
//	 6    空 0    0    空表示这一维是自由边界  其它0处，可以是0或非零值 但必须有值
//   5     0   空 0    
//	 4   空 空 0
//   3     0    0  空
//   2   空 0   空
//   1    0   空 空
{
	int IA,IB,IH,LK;

	for( int I=1; I <= M14; I++  )
	{
		IA=2*(I-1);
		IB=3*(I-1);
		if( IP[IA] > N ) goto L75;
		IH=3*(IP[IA]-1);
		LK=IP[IA+1];
		switch( LK )
		{
		case 1:	goto L35;
		case 2:	goto L40;
		case 3:	goto L45;
		case 4:	goto L50;
		case 5:	goto L55;
		case 6:	goto L60;
		case 7:	goto L65;
		default: break;
		}
L35:		if( VP[IB] == 0.0 ) goto L75;
		if( KG == 0 ) goto L36;
		W[IH]=VP[IB];
		goto L75;
L36:		if( KG1 == 0 ) goto L136;
		G[IH]=0.0;
L136:	W[IH]=0.0;
		goto L75;
L40:		if( VP[IB+1] == 0.0 ) goto L75;
		if( KG == 0 ) goto L41;
		W[IH+1]=VP[IB+1];
		goto L75;
L41:		if( KG1 == 0 ) goto L141;                                            
		G[IH+1]=0.0;
L141:	W[IH+1]=0.0;
		goto L75;
L45:		if( VP[IB] == 0.0 ) goto L47;
		if( KG == 0) goto L46;                                              
		W[IH]=VP[IB];                                                 
		goto L47;
L46:		if( KG1 == 0 ) goto L146;
		G[IH]=0.0;
L146:	W[IH]=0.0;
L47:		if( VP[IB+1] == 0.0 ) goto L75;
		if( KG == 0) goto L48;
		W[IH+1]=VP[IB+1];                                                  
		goto L75;
L48:		if( KG1 == 0 ) goto L148;                                            
		G[IH+1]=0.0;
L148:	W[IH+1]=0.0;
		goto L75;
L50:		if( VP[IB+2] == 0.0 ) goto L75;
		if( KG == 0 ) goto L51;                                             
		W[IH+2]=VP[IB+2];                                                  
		goto L75;                                                          
L51:		if( KG1 == 0 ) goto L151;                                            
		G[IH+2]=0.0;
L151:	W[IH+2]=0.0;
		goto L75;
L55:		if( VP[IB] == 0.0 ) goto L57;
		if( KG == 0) goto L56;                                              
		W[IH]=VP[IB];                                                 
		goto L57;
L56:		if(KG1 == 0) goto L156;                                            
		G[IH]=0.0;
L156:	W[IH]=0.0;
L57:		if( VP[IB+2] == 0.0 ) goto L75;
		if( KG == 0 ) goto L58;                                              
		W[IH+2]=VP[IB+2];
		goto L75;
L58:		if( KG1 == 0 ) goto L158;                                            
		G[IH+2]=0.0;
L158:	W[IH+2]=0.0;
		goto L75;
L60:		if( VP[IB+1] == 0.0 ) goto L62;
		if( KG == 0 ) goto L61;                                              
		W[IH+1]=VP[IB+1];                                                  
		goto L62;
L61:		if( KG1 == 0 ) goto L161;                                            
		G[IH+1]=0.0;
L161:	W[IH+1]=0.0;
L62:		if( VP[IB+2] == 0.0) goto L75;
		if( KG == 0 ) goto L63;                                              
		W[IH+2]=VP[IB+2];
		goto L75;
L63:		if( KG1 == 0) goto L163;                                            
		G[IH+2]=0.0;
L163:	W[IH+2]=0.0;
		goto L75;
L65:		if( VP[IB] == 0.0 ) goto L67;
		if( KG == 0 ) goto L66;                                              
		W[IH]=VP[IB];
		goto L67;
L66:		if( KG1 == 0) goto L166;
		G[IH]=0.0;
L166:	W[IH]=0.0;
L67:		if( VP[IB+1] == 0.0 ) goto L69;
		if( KG == 0) goto L68;
		W[IH+1]=VP[IB+1];
		goto L69;
L68:		if( KG1 == 0) goto L168;                                            
		G[IH+1]=0.0;
L168:	W[IH+1]=0.0;
L69:		if( VP[IB+2] == 0.0 ) goto L75;
		if( KG == 0 ) goto L70;                                              
		W[IH+2]=VP[IB+2];                                                  
		goto L75;
L70:		if( KG1 == 0 ) goto L170;                                            
		G[IH+2]=0.0;
L170:	W[IH+2]=0.0;
L75:	;
	}
}

//---------------------------------------------------------------------------
void SolveEqu::mabvm( int N, int N1, int *IZ, int *IG, double *AK, double *U, double *V )
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
//     CALLED BY SOL.																		*     
//     PURPOSE:																					*     
//     MATRIX AK BY VECTOR U MULTIPLICATION.						*     
//     N_TOTAL OF NODES.																*     
//     N1_2*N OR 3*N.																			*     
//     IZ_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.	*     
//     IG_AUXILIARY INFORMATION FOR STIFFNESS MATRIX.	*     
//     AK_STIFFNESS MATRIX.															*     
//     U_VECTOR.																					*     
//     V_RESULT STORAGE.																*     
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
{
	int IH,IH1,IH2,IK,KH,IC,IA;

	for( int I=0; I < N1; I++  )
	{
		V[I]=0.0;
	}

	IH=3;                                                              
	IH1=6;                                                             
	IH2=9;
	for( int I=1; I <= N; I++  )
	{
		IK=IH2*IZ[I-1]+IH1*(I-1);                                            
		KH=IH*(I-1);
		for( int J=1; J <= IH; J++ )
		{
			for( int K=1; K <= IH; K++ )
			{
				if(J == 3) goto L25;                                               
				if(K == 3) goto L25;                                               
				V[KH+J-1]=V[KH+J-1]+AK[IK+J+K-2]*U[KH+K-1];                              
				goto L30;
L25:				V[KH+J-1]=V[KH+J-1]+AK[IK+J+K-1]*U[KH+K-1];
L30:	;
			}
		}
		if(I <= 1) goto L60;
		IC=IZ[I-1]-IZ[I-2];
		for( int L=1; L <= IC; L++ )
		{
			IK=IH1*(I-1)+IH2*(IZ[I-2]+L-1);
			IA=IH*(IG[IZ[I-2]+L-1]-1);
			for( int J=1; J <= IH; J++ )
			{
				for( int K=1; K <= IH; K++ )
				{
					V[IA+J-1]=V[IA+J-1]+AK[IK+IH*(K-1)+J-1]*U[KH+K-1];                         
					V[KH+J-1]=V[KH+J-1]+AK[IK+IH*(J-1)+K-1]*U[IA+K-1];                        
				}
			}
		}
L60:		;
	}
}

//===========================================================================
