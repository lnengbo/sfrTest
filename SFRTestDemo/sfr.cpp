#include "sfr.h"
#include "opencv2/opencv.hpp"
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <QDebug>

void de_Gamma(cv::Mat &Src, double gamma)
{
	if (Src.channels() != 1) { return; }

	for (int i = 0; i < Src.rows; ++i)
	{
		uchar *SrcP = Src.ptr(i);
		for (int j = 0; j < Src.cols; ++j)
		{
			SrcP[j] = 255 * (pow((double)SrcP[j] / 255, 1 / gamma));
		}
	}
}

void SLR(std::vector<double> &Cen_Shifts, std::vector<double> &y_shifts, double *a, double *b)
{
	//a -> intercept, b->slope of slanted edge

	int y_size = y_shifts.size();
    //用简单线性回归求解方程，得到斜率和截距
	double xsquare = 0, xavg = 0, yavg = 0;
	int i;
	for (i = 0; i < y_size; ++i)
	{
		yavg += y_shifts[i];
		xavg += Cen_Shifts[i];
	}
	xavg /= (double)y_size;
	yavg /= (double)y_size;

    //  简单线性回归 Cen_Shifts//x
	for (i = 0; i < y_size; ++i)
	{
		double temp = (Cen_Shifts[i] - xavg);
        *b += temp * (y_shifts[i]); //(y_shifts[i] - yavg)
		xsquare += temp * temp;
	}

    *b /= xsquare;//斜率
    *a = yavg - (*b)*xavg;//截距
    qDebug() << "rate" << *b;
}

std::vector<double> CentroidFind(cv::Mat &Src,std::vector<double> &y_shifts, double *CCoffset)
{
	std::vector<double> Cen_Shifts(Src.rows);
    int i, j, height = Src.rows, width = Src.cols;
	cv::Mat tempSrc(Src.size(), CV_8UC1);
    cv::Mat copySrc(tempSrc.size(),CV_8UC1);
    cv::medianBlur(Src,Src,5);
    Src.copyTo(copySrc);
    double molecule = 0, denominator = 0, temp = 0;
    double minVal = 0.0;
    double maxVal = 0.0;
    cv::Point  minLoc;
    cv::Point  maxLoc;
    cv::minMaxLoc(Src,&minVal,&maxVal,&minLoc,&maxLoc);
//    if(minVal / maxVal == 0)
//        return {};
    //对模板感兴趣区域图像进行双边滤波，以确保找到更准确的倾斜边缘
    Src.copyTo(tempSrc);
    //cv::imshow("tempSrc",tempSrc);
    // 计算ROI中每一行的质心
    // 形心模板：总力矩/总量  差分
    for (int i = 0; i < height; ++i)
	{
        molecule = 0.0;
        denominator =0.0;
        temp = 0.0;
		uchar *tempSrcPtr = tempSrc.ptr<uchar>(i);
        for (int j = 0; j < width -1; j++)
		{
            temp = (double)tempSrcPtr[j + 1] - (double)tempSrcPtr[j] ;
			molecule += temp * j;
			denominator += temp;
		}
        Cen_Shifts[i] = (double)molecule / denominator ;
        copySrc.at<uchar>(i,Cen_Shifts[i]) = 255;
	}

    //消除远离斜边的噪声(+/- 10 pixels)
	tempSrc = Src.clone();
    //cv::GaussianBlur(Src, Src, cv::Size(3, 3), 0);

    //替换去除周边像素
    for (i = 0; i < height; ++i)
    {
        uchar *SrcPtr = Src.ptr<uchar>(i);
        uchar *tempPtr = tempSrc.ptr<uchar>(i);
        for (j = int(Cen_Shifts[i]) - 10; j < int(Cen_Shifts[i]) + 10; ++j)
        {
            SrcPtr[j] = tempPtr[j];
        }
    }
    // 检查图像中的边缘是否过于靠近图像的角点
	if (Cen_Shifts[0] < 2.0 || width - Cen_Shifts[0] < 2.0)
	{
		std::cerr << "The edge in ROI is too close to the image corners" << std::endl;
        return std::vector<double>(0,0);
	}

	if (Cen_Shifts[height - 1] < 2 || width - Cen_Shifts[height - 1] < 2)
	{
		std::cerr << "The edge in ROI is too close to the image corners" << std::endl;
        return std::vector<double>(0,0);
	}

    double half_ySize = height / 2;
    double CC = Cen_Shifts[half_ySize]; //像心质心

    for (i = 0; i < height; ++i)
    {
        Cen_Shifts[i] -= CC; //计算质心和图像中心质心之间的偏移
        y_shifts[i] = (double)i - half_ySize ; //计算图像中心高度和每行之间的偏移

    }
    cv::resize(copySrc,copySrc,cv::Size(300,300));
    //cv::imshow(QString("copySrc%1").arg(Cen_Shifts[0]).toLocal8Bit().data(),copySrc);
    *CCoffset = CC;
	return Cen_Shifts;
}


std::vector<double> OverSampling(cv::Mat &Src, double slope, double CCoffset, int height, int width, int *SamplingLen,std::vector<double> & esfVec)
{
	std::vector<double> RowShifts(height, 0);

	int i, j, k;
    int halfY = height /2;


    //计算每行的像素偏移，使每行的质心尽可能接近 基于拟合结果更新shifts 转成ESF线段
    //slope斜率
    for (i = 0; i < height; ++i) {
        RowShifts[i] = ((double)(i - halfY) / slope) + CCoffset;
        //std::cout<< "ROW" << RowShifts[i] << std::endl;
    }


    //DataMap->记录移位后像素的索引

    //数据->记录原始像素值
	std::vector<double> DataMap(height*width, 0);
	std::vector<double> Datas(height*width, 0);

    //数据归一化
	for (i = 0, k = 0; i < height; ++i)
	{
        int baseindex = width * i;
		for (j = 0; j < width; ++j)
		{
            DataMap[baseindex + j] = j - RowShifts[i]; ////计算每个点到刀口的距离
            Datas[baseindex + j] = int(Src.at<uchar>(i, j) - 0) / (1 - 0);
        }
	}

	std::vector<double> SamplingBar(*SamplingLen, 0);
	std::vector<int> MappingCount(*SamplingLen, 0);
    //开始将原始数据映射到4x采样数据，并记录4x采样数据中每个像素的计数
	for (i = 0; i < height*width; ++i)
	{
        int MappingIndex = static_cast<int>(4 * DataMap[i]);//向下取整
		if (MappingIndex >= 0 && MappingIndex < *SamplingLen)
		{
			SamplingBar[MappingIndex] = SamplingBar[MappingIndex] + Datas[i];
			MappingCount[MappingIndex]++;
		}
	}
    int nzeros = 0;
    //取4x采样数据中的像素值平均，如果像素中的像素值为零，则复制接近像素的值
    for (i = 0; i < *SamplingLen; ++i) {
		j = 0;
		k = 1;
		if (MappingCount[i] == 0) {
            nzeros++;
			if (i == 0) {
                while (!j)//当j == 0 时，表示此处信号为0
				{
                    if (MappingCount[i + k] != 0)//第一行元素
					{
						SamplingBar[i] = SamplingBar[i + k] / ((double)MappingCount[i + k]);
                        j = 1;//充当Flag;
					}
                    else ++k;//找到不为0
				}
			}
			else {
				while (!j && ((i - k) >= 0))
				{
                    //j==0&&i-k>=0  j==0说明 counts[i]是0  i-k>0说明  k在i前面，找前面不为零的数值赋给AveEdge[i]
					if (MappingCount[i - k] != 0)
					{
						SamplingBar[i] = SamplingBar[i - k];   /* Don't divide by counts since it already happened in previous iteration */
						j = 1;
					}
					else ++k;
				}
				if ((i - k) < 0)
				{
                    ////k>i,其实联合上面那段，就是：此处的counts[i]累加次数为零，所以AveEdge[i]也就是0，所以要找到附近一个不为零的值赋给AveEdge[i]
					k = 1;
					while (!j)
					{
						if (MappingCount[i + k] != 0)
						{
							SamplingBar[i] = SamplingBar[i + k] / ((double)MappingCount[i + k]);
							j = 1;
						}
						else ++k;
					}
				}
			}
		}
		else
			SamplingBar[i] = (SamplingBar[i]) / ((double)MappingCount[i]);
        ////如果此处不为零，直接就求个平均值
	}

    //减少采样数据的长度
    //由于边缘数据是唯一重要的数据，所以我们将靠近采样数据两侧的数据截断
    //对结果的贡献很小
    //截断小于原始长度10%且大于原始长度90%的数据
	int originalSamplingLen = *SamplingLen;
    *SamplingLen = *SamplingLen * 0.9;

    if (nzeros > 0) //提示信息
        std::cout << "\nWARNING: "<< nzeros << "Zero counts found during projection binning.\n" << std::endl;

     std::vector<double> deSampling(*SamplingLen,0);
    for(int i =0;i < SamplingBar.size();i++)
        esfVec.push_back(SamplingBar[i]);
    //导数采样数据（ESF）得到线扩展函数 差分运算
    for (i = originalSamplingLen * 0.1, j = 1; i < originalSamplingLen, j < *SamplingLen; ++i,++j)
    {
        deSampling[j] = SamplingBar[i + 1] - SamplingBar[i];
        //std::cout << "ESFSample" << deSampling[j];
    }

	return deSampling;
}


std::vector<double> HammingWindows(std::vector<double> &deSampling, int SamplingLen)
{
	int i, j;
	std::vector<double> tempData(SamplingLen);

    //我们希望将峰值数据转移到线扩散函数数据的中心

    //因为我们稍后会做汉明窗口，这将使重要数据远离过滤

    //如果有两个峰值，我们使用两个变量来记录峰值数据的位置
	int L_location = -1, R_location = -1;
    long begin = 0;
	double SamplingMax = 0;
	for (i = 0; i < SamplingLen; ++i)
	{
		if (fabs(deSampling[i]) > fabs(SamplingMax))
			SamplingMax = deSampling[i];
	}

	for (i = 0; i < SamplingLen; ++i)
	{
		if (deSampling[i] == SamplingMax)
		{
            if (L_location < 0) {
                L_location = i;
            }
			R_location = i;
		}
	}

    //移位量
	int PeakOffset = (R_location + L_location) / 2 - SamplingLen / 2;
	if (PeakOffset)
	{
		for (i = 0; i < SamplingLen; ++i)
		{
			int newIndex = i - PeakOffset;
			if (newIndex >= 0 && newIndex < SamplingLen) { tempData[newIndex] = deSampling[i]; }
		}
	}
	else
	{
		for (i = 0; i < SamplingLen; ++i) { tempData[i] = deSampling[i]; }
	}
    int x = begin;
    int y = -SamplingLen/2;
    //汉明窗过滤
    for (; x< SamplingLen; ++x,++y)
	{
        tempData[x] = tempData[x] * (0.54 + 0.46*cos(2 * M_PI*(double)y / (SamplingLen / 2)));
	}

	return tempData;
}

void DFT(std::vector<double> &data, int size)
{
	int i, j;
	std::complex<double> *arr = new std::complex<double>[size];
	for (i = 0; i < size; ++i)
	{
		arr[i] = std::complex<double>(data[i], 0);
	}

	for (i = 0; i < size / 2.0; ++i)
	{
		std::complex<double> temp = 0;
		for (j = 0; j < size; ++j)
		{
            double w = 2 * M_PI * i * j / size;
			std::complex<double> deg(cos(w), -sin(w));
			temp += arr[j] * deg;
		}
		data[i] = sqrt(temp.real()*temp.real() + temp.imag()*temp.imag());
	}
}

void ReduceRows(double slope, int *ImgHeight)
{
	double tempSlope = fabs(slope);
    int numcycles = (*ImgHeight) * tempSlope;
    double cyclelimit = 5.0;
    if(tempSlope < cyclelimit/(double)(*ImgHeight))
    {
        std::cout << "Edget angle" <<atan(tempSlope) * 180 /M_PI <<" is so shallow it needs\n" << std::endl;
        return;
    }

    if ( numcycles/tempSlope <= *ImgHeight)
    {
        *ImgHeight = (numcycles / tempSlope);
    }
}

int SFRCalculation(cv::Mat &ROI, double gamma,std::vector<double> &resVec,std::vector<double> &esfVec,std::vector<double> &lsfVec)
{
	if (ROI.empty())
	{
		std::cerr << "Open the ROI image error" << std::endl;
		return 0;
	}
	int height = ROI.rows, width = ROI.cols;
    //进行gamma解码以消除摄像机设备编码的gamma
	de_Gamma(ROI, gamma);
	int i, j;

    double slope = 0, intercept = 0;
	
    //中心质心偏移
	double CCoffset = 0;
	std::vector<double> y_shifts(height);

	////Calculate the shifts between Centroids and Centroid of image center	
    /// 计算中心
	std::vector<double> Cen_Shifts = CentroidFind(ROI, y_shifts, &CCoffset);
    qDebug() << "Cen" << Cen_Shifts[10] << CCoffset;

    if (Cen_Shifts.empty()) {
        return 0;
    }

    //斜边拟合的简单线性回归
	SLR(Cen_Shifts, y_shifts, &intercept, &slope);

    //将数据行数截断为最大斜率周期，该周期将具有整数个全相位旋转
	ReduceRows(slope, &height);

    //将cOffset更新为原始图像中点和我们计算的参考中点之间的偏移
	CCoffset = CCoffset + 0.5 + intercept - width / 2;

    CCoffset += 0.155;

    //将原始图像的像素值映射为长度为原始图像宽度4倍的采样数据

    //此步骤用于集中原始像素值的变化量
	int SamplingLen = width * 4;
    std::vector<double> OverSamplingData = OverSampling(ROI, slope, CCoffset, height, width, &SamplingLen,esfVec);

    for(int i = 0; i < OverSamplingData.size();i++)
        lsfVec.push_back(OverSamplingData[i]);
    //利用汉明窗对数据两侧的纹波信号进行滤波
     OverSamplingData = HammingWindows(OverSamplingData, SamplingLen);

    //四变换
	DFT(OverSamplingData, SamplingLen);
	width = int(SamplingLen / 4);
	double maxData = 0;
	for (i = 0; i < SamplingLen; ++i)
	{
        if (OverSamplingData[i] != 0)
		{

			maxData = OverSamplingData[i];
            std::cout << "maxData " << maxData << std::endl;
            break;
		}
	}

    for (int i = 0; i < SamplingLen; ++i)
    {
        OverSamplingData[i] /= maxData;
        //std::cout << "Oversample" <<OverSamplingData[i] << std::endl;
    }

	for (i = 0; i <= width; ++i)
	{
		double frequency = (double)i / width;
        resVec.push_back(abs(OverSamplingData[i]));
        if(i == width /2)
            qDebug() << frequency << "," << OverSamplingData[i];
	}

	return 1;
}
