#include "local_cartesian.h"

namespace proj
{
	LocalCartesian::LocalCartesian()
	{
	}

	LocalCartesian::~LocalCartesian()
	{
	}

	void LocalCartesian::set_home(double lon, double lat, double alt)
	{
		home_lon_ = lon;
		home_lat_ = lat;
		home_alt_ = alt;
	}

	bool LocalCartesian::LLATOENU(const double lon,const double lat,const double alt, double &x, double &y, double &z)
	{
		if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || alt < -3000.0 || alt > 10000.0)
			return false;

		Eigen::Vector3d lla{lon, lat, alt};
		Eigen::Vector3d home_lla{home_lon_, home_lat_, home_alt_};
		Eigen::Vector3d ecef_xyz{0, 0, 0};
		Eigen::Vector3d base_xyz{0, 0, 0};
		Eigen::Vector3d enu{0, 0, 0};

        //旋转矩阵 ECEF坐标系 先绕Z轴旋转 90+Lon , 再绕X轴 旋转 90-Lat 最后 平移可到达 ENU坐标系
		//坐标变换 ： ENU 坐标系下坐标  = 旋转矩阵逆矩阵 x ECEF坐标系下坐标
	    Eigen::AngleAxisd Z_V(-(M_PI_2 + home_lla(0) * DEG2RAD), Eigen::Vector3d(0, 0, 1));
		Eigen::Matrix3d Z_R = Z_V.matrix();

        Eigen::AngleAxisd X_V(home_lla(1) * DEG2RAD - M_PI_2, Eigen::Vector3d(1, 0, 0));
		Eigen::Matrix3d X_R = X_V.matrix();		

        Eigen::Matrix<double, 3, 3> Matrix_R = X_R * Z_R; //逆矩阵 右乘组合

        //先将地理坐标系下坐标 统一转换至ECEF坐标系坐标
		/*************** HOME(LLA) TO ECEF *******************/
		LLATOECEF(home_lla,base_xyz);

        /****************** LLA TO ECEF *********************/
		LLATOECEF(lla,ecef_xyz);

       //后将ECEF坐标系 坐标变换至 ENU坐标系坐标
       /****************** ECEF TO ENU *********************/
	   //先平移后 通过旋转矩阵逆矩阵 进行坐标转换
       enu  =  Matrix_R * (ecef_xyz -  base_xyz);

	   x = enu(0);
	   y = enu(1);
	   z = enu(2);

	   if (std::isnan(x)|| std::isinf(x) || std::isnan(y)|| std::isinf(y) || std::isnan(z)|| std::isinf(z))
	   {
		x = 0;
		y = 0;
		z = 0;
		return false;
	   }

	   return true;
	}

	bool LocalCartesian::ENUTOLLA(const double x,const double y,const double z, double &lon, double &lat, double &alt)
	{
		
		Eigen::Vector3d enu{x, y, z};
		Eigen::Vector3d home_lla{home_lon_, home_lat_, home_alt_};
		Eigen::Vector3d ecef_xyz{0, 0, 0};
		Eigen::Vector3d base_xyz{0, 0, 0};
        
        //旋转矩阵 ECEF坐标系 先绕Z轴旋转 90+Lon , 再绕X轴 旋转 90-Lat 最后 平移可到达 ENU
		//坐标变换 ：ECEF坐标系下坐标  = 旋转矩阵 x ENU坐标系下坐标 
        Eigen::AngleAxisd Z_V(M_PI_2 + home_lla(0) * DEG2RAD, Eigen::Vector3d(0, 0, 1));
		Eigen::Matrix3d Z_R = Z_V.matrix();

        Eigen::AngleAxisd X_V(M_PI_2 - home_lla(1) * DEG2RAD, Eigen::Vector3d(1, 0, 0));
		Eigen::Matrix3d X_R = X_V.matrix();		

        Eigen::Matrix<double, 3, 3> Matrix_R = Z_R * X_R; //左乘

        //HOME点转换至ECEF坐标
		/****************** HOME(LLA) TO ECEF *********************/
		LLATOECEF(home_lla,base_xyz);
        
		//ENU坐标系坐标 先进行旋转 后 平移 得到 ECEF坐标系坐标
	    /****************** ENU TO ECEF *********************/
        ecef_xyz = Matrix_R * enu + base_xyz;

        //ECEF坐标系转换至经纬度高地理坐标系有两种方法 迭代求解法(iterative method)和 闭式解(a closed form solution)
		//哪个更快取决于编程环境、计算架构以及需要多少精度。
		//因为两种方法的计算量大致相同，但第一种方法必须至少迭代一次，所以实际上可能会慢一些。
		//这里采用闭式解
		/****************** ECEF TO LLA *********************/
		double P = sqrtf(ecef_xyz[0] * ecef_xyz[0] + ecef_xyz[1] * ecef_xyz[1]);
        double U = atan2(ecef_xyz[2] * A_Earth , B_Earth * P);

        lon = atan2(ecef_xyz[1],ecef_xyz[0]) * RAD2DEG;
		lat = atan2(ecef_xyz[2]+ e2_2 * B_Earth * pow(sin(U),3),P - A_Earth * e1_2 * pow(cos(U),3)) * RAD2DEG;
		double N = A_Earth / (sqrt(1 - e1_2 * pow(sin(lat * DEG2RAD),2)));

        double cos_lat = limitEps(cos(lat * DEG2RAD),1e-5);
		alt = P / cos_lat - N;

		if (std::isnan(lon)|| std::isnan(lat)|| std::isnan(alt) ||std::isinf(lon)||std:: isinf(lat)|| std::isinf(alt))
		{
			lon = 0;
			lat = 0;
			alt = 0;
			return false;
		}

		return true;
	}

	void LocalCartesian::LLATOECEF(const Eigen::Vector3d lla,Eigen::Vector3d& ecef)
	{
		double sin_lon = sin((lla(0) * DEG2RAD));
		double cos_lon = cos((lla(0) * DEG2RAD));
		double sin_lat = sin((lla(1) * DEG2RAD));
		double cos_lat = cos((lla(1) * DEG2RAD));

		double N = A_Earth / (sqrtf(1 - e1_2 * sin_lat * sin_lat));
		ecef(0) = (N + lla(2)) * cos_lat * cos_lon;
		ecef(1) = (N + lla(2)) * cos_lat * sin_lon;
		ecef(2) = (N * (1 - e1_2) + lla(2))  * sin_lat;
	}
} // namespace proj