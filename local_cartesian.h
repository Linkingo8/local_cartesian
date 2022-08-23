/**
 * @file local_cartesian.h
 * @author linheng
 * @brief 经纬度与ENU坐标投影变换
 * @version 1.0
 * @date 2022-08-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef LOCAL_CARTESIAN_H
#define LOCAL_CARTESIAN_H

#include <cmath>
#include <iostream>
#include <array>
#include <float.h>

#include <eigen3/Eigen/Dense>


# define RAD2DEG	57.295779513082
# define DEG2RAD  0.0174532925199

namespace proj
{
  class LocalCartesian
  {
  private:
    // home点的经度，纬度，海拔高度
    double home_lon_, home_lat_, home_alt_;

    // 坐标系参数，WGS84
    const double A_Earth = 6378137.0;           // 长半径
    const double B_Earth = 6356752.3142;        // 短半径
    const double F_Earth = 1.0 / 298.257223563; // 扁率

    const double e1_2 = 0.00669437999013;  // 第一偏心率平方
    const double e2_2 = 0.00673949674227;  // 第二偏心率平方

    /**
     * @brief 经纬度转ECEF
     * 
     * @param lon 经度
     * @param lat 纬度
     * @param alt 高度
     * @param x 
     * @param y 
     * @param z 
     * @return true 
     * @return false 
     */
    void LLATOECEF(const Eigen::Vector3d lla,Eigen::Vector3d& ecef);
    
    //限制最小值
    double limitEps(double val, const double eps)
    {
        if(val >= 0 && val < eps)
            val = eps;
        else if(val <= 0 && val > -eps)
            val = -eps;
        return val;
    }
  public:
    LocalCartesian();
    ~LocalCartesian();

    /**
     * @brief 设置home点经纬度，高度
     * 
     * @param lon home点经度
     * @param lat home点纬度
     * @param alt home点高度
     */
    void set_home(double lon, double lat, double alt);

    /**
     * @brief 经纬度转local(ENU)
     * 
     * @param lon 经度
     * @param lat 纬度
     * @param alt 高度
     * @param x 
     * @param y 
     * @param z 
     * @return true 
     * @return false 
     */
    bool LLATOENU(const double lon,const double lat,const double alt, double &x, double &y, double &z);

    /**
     * @brief local(ENU)转经纬度
     *        
     * @param x 
     * @param y 
     * @param z 
     * @param lon 经度
     * @param lat 纬度
     * @param alt 高度
     * @return true 
     * @return false 
     */
    bool ENUTOLLA(const double x,const double y,const double z, double &lon, double &lat, double &alt);
  };
} // namespace proj

#endif