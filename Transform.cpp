#include"Transform.h"
#include "MyMath.h"
#include<cmath>

Transform::Transform() {
	width = height = 0.0;
}

Transform::~Transform() {
	width = height = 0;
}


//返回一个右手坐标系的透视投影矩阵，类似D3DXMatrixPerspectiveFovRH,参考https://msdn.microsoft.com/en-us/library/windows/desktop/bb205351(v=vs.85).aspx
//fovy：观察时y轴方向的角度（弧度）
//aspect：纵横比，视线空间的宽度除以高度
//zn:近裁剪面的Z值
//zf：远裁剪面的Z值
void Transform::set_projection_matrix(double fovy, double aspect, double zn ,double zf) {
	/*double theCot = 1.0f / tan(fovy / 2.0f);
	this->projection_matrix.m[0][0] = theCot / aspect;
	this->projection_matrix.m[1][1] = theCot;
	this->projection_matrix.m[2][2] = zf / (zf - zn);
	this->projection_matrix.m[2][3] = -1;
	this->projection_matrix.m[3][2] = zn*zf / (zf - zn);*/

	double yScale = 1.0f / tan(fovy / 2.0f);
	double xScale = yScale / aspect;
	this->projection_matrix.m[0][0] = xScale;
	this->projection_matrix.m[1][1] = yScale;
	this->projection_matrix.m[2][2] = zf / (zf - zn);
	this->projection_matrix.m[2][3] = -1;
	this->projection_matrix.m[3][2] = zn * zf / (zf - zn);



}


//视角变换矩阵，类似于LookAt，参考https://msdn.microsoft.com/en-us/library/windows/desktop/bb281711(v=vs.85).aspx
//cameraposition : 摄像机的位置
//cameraTarget：摄像机朝向的方向
//up：摄像机向上的方向

void Transform::set_view_matrix(vect* cameraPosition, vect* cameraTarget, vect* up) {
	vect zaxis, xaxis, yaxis, subvect;
	vect_sub(&subvect, cameraPosition, cameraTarget);
	vect_normalize(&zaxis, &subvect);

	vect_crossmul(&subvect, up, &zaxis);
	vect_normalize(&xaxis, &subvect);
	vect_crossmul(&yaxis, &zaxis, &xaxis);

	this->view_matrix.m[0][0] = xaxis.x;
	this->view_matrix.m[0][1] = yaxis.x;
	this->view_matrix.m[0][2] = zaxis.x;
	this->view_matrix.m[0][3] = 0;

	this->view_matrix.m[1][0] = xaxis.y;
	this->view_matrix.m[1][1] = yaxis.y;
	this->view_matrix.m[1][2] = zaxis.y;
	this->view_matrix.m[1][3] = 0;

	this->view_matrix.m[2][0] = xaxis.z;
	this->view_matrix.m[2][1] = yaxis.z;
	this->view_matrix.m[2][2] = zaxis.z;
	this->view_matrix.m[2][3] = 0;

	this->view_matrix.m[3][0] = -vect_dotmul(&xaxis, cameraPosition);
	this->view_matrix.m[3][1] = -vect_dotmul(&yaxis, cameraPosition);
	this->view_matrix.m[3][2] = -vect_dotmul(&zaxis, cameraPosition);
	this->view_matrix.m[3][3] = 1;
	//
	//matr T_view;
	//matr_identity(&T_view);
	//matr R_view;
	//matr_identity(&R_view);
	//vect temp,g;
	//vect_sub(&temp, cameraPosition, cameraTarget);
	//vect_normalize(&g, &temp);

	//T_view.m[0][3] = -cameraPosition->x;
	//T_view.m[1][3] = -cameraPosition->y;
	//T_view.m[2][3] = -cameraPosition->z;
	//

	//vect gt;
	//vect_crossmul(&temp, &g, up);
	////vect_normalize(&gt, &temp);
	//R_view.m[0][0] = gt.x;
	//R_view.m[0][1] = gt.y;
	//R_view.m[0][2] = gt.z;
	//R_view.m[1][0] = up->x;
	//R_view.m[1][1] = up->y;
	//R_view.m[1][2] = up->z;
	//R_view.m[2][0] = g.x;
	//R_view.m[2][1] = g.y;
	//R_view.m[2][2] = g.z;

	//matrix_mul(&this->view_matrix, &R_view, &T_view);


}

//得到最终的转换矩阵
void Transform::set_transform_matrix() {
	matr tempmatr;
	matrix_mul(&tempmatr, &this->world_matrix, &this->view_matrix);
	matrix_mul(&this->transform_matrix, &tempmatr, &this->projection_matrix);
}

void getScreenPos(vect* out, vect* v_in, double width, double height) {
	double rate = 1.0f / v_in->w;
	out->x = (int)((v_in->x * rate + 1.0f) * width * 0.5f + 0.5f);
	out->y = (int)((1.0f - v_in->y * rate) * height * 0.5f + 0.5f);
	out->z = v_in->z * rate;
	out->w = 1.0f;
}
