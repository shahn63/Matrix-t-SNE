#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gradient(NumericMatrix mat1, NumericMatrix mat2, NumericMatrix mats, int n1, int n2) {

  // 결과를 저장할 행렬
	NumericMatrix value1(n1,2);
	NumericMatrix mat3(n1,2);
	NumericMatrix grads1(n1,n2), grads2(n1,n2), grads(n1,n2+n2);

  // 벡터의 원소를 제곱하고 더하기
	for (int jj = 0; jj < n2; jj++){
		mat3(_,0)=mat1(_,jj);
		mat3(_,1)=mat2(_,jj);
		for (int i = 0; i < n1; i++){
			for (int j=0; j < n1; j++){
       				value1(j, 0) = -mat3(j, 0) + mat3(i, 0);
        			value1(j, 1) = -mat3(j, 1) + mat3(i, 1);
        			value1(j, 0) *= mats(j, i);
		        	value1(j, 1) *= mats(j, i);
			}
			for (int j=0; j < n1; j++){
				grads1(i,jj) += value1(j,0);
				grads2(i,jj) += value1(j,1);
			}
		}
	}

  	for (int i=0; i< n2; i++){
		grads(_,i)=grads1(_,i);
	}
	for (int i=n2; i<n2+n2; i++){
		grads(_,i)=grads2(_,i-n2);
	}	


	return grads;
}
