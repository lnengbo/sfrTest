#include "mainwindow.h"
#include <QApplication>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <QDebug>
#include <QTime>
#include "qcustomplot.h"
QCustomPlot *custplot;
QCustomPlot *custplotesf;
QCustomPlot *custplotlsf;
cv::Mat src, gray_src;
char *output_image =(char*) "Tomoshi Result";		//输出窗口
int num_corners = 1;		//角点个数
const int max_corners = 50;		//定义最大的角点个数

/*---------------------------------------------------------------------*/
#include "sfr.h"
void ROImouseEvent(int event, int x, int y, int flag, void *params)
{
    cv::Point *ptr = (cv::Point*)params;

    if (event == cv::EVENT_LBUTTONDOWN && ptr[0].x == -1 && ptr[0].y == -1)
    {
        ptr[0].x = x;
        ptr[0].y = y;
    }

    if (flag == cv::EVENT_FLAG_LBUTTON)
    {
        ptr[1].x = x;
        ptr[1].y = y;
    }

    if (event == cv::EVENT_LBUTTONUP && ptr[2].x == -1 && ptr[2].y == -1)
    {
        ptr[2].x = x;
        ptr[2].y = y;
    }

}

void ROISelection(cv::Mat &img,cv::Mat &roi)
{
    cv::Point *Corners = new cv::Point[3];
    Corners[0].x = Corners[0].y = -1;
    Corners[1].x = Corners[1].y = -1;
    Corners[2].x = Corners[2].y = -1;

    cv::namedWindow("ROI select(Press Esc to close window)",cv::WINDOW_NORMAL);
    cv::imshow("ROI select(Press Esc to close window)", img);

    bool downFlag = false, upFlag = false;
    while (cv::waitKey(1) != 27)
    {
        cv::setMouseCallback("ROI select(Press Esc to close window)", ROImouseEvent, Corners);

        if (Corners[0].x != -1 && Corners[0].y != -1) { downFlag = true; }
        if (Corners[2].x != -1 && Corners[2].y != -1) { upFlag = true; }

        if (downFlag && !upFlag && Corners[1].x != -1)
        {
            cv::Mat LocalImage = img.clone();
            cv::rectangle(LocalImage, Corners[0], Corners[1], cv::Scalar(255, 255, 255), 2);
            cv::imshow("ROI select(Press Esc to close window)", LocalImage);
        }

        if (downFlag && upFlag)
        {

            cv::Rect ROIBox;
            ROIBox.width = abs(Corners[0].x - Corners[2].x);
            ROIBox.height = abs(Corners[0].y - Corners[2].y);

            if (ROIBox.width < 5 && ROIBox.height < 5)
            {
                std::cerr << "ROI size too small, please re-crop the ROI" << std::endl;
            }


            ROIBox.x = Corners[0].x < Corners[1].x ? Corners[0].x : Corners[1].x;
            ROIBox.y = Corners[0].y < Corners[1].y ? Corners[0].y : Corners[1].y;

            roi = img(ROIBox);
            downFlag = upFlag = false;

            Corners[0].x = Corners[0].y = -1;
            Corners[1].x = Corners[1].y = -1;
            Corners[2].x = Corners[2].y = -1;

        }

    }
    cv::destroyWindow("ROI select(Press Esc to close window)");

    delete[] Corners;
}
void Tomoshi_demo(int, void*) {
    if (num_corners < 2)
        num_corners = 1;		//避免角点太少

    std::vector<cv::Point2f> corners;
    double qualityLevel = 0.1;		//最大，最小特征值的乘法因子。定义可接受图像角点的最小质量因子。
    double min_distance = 20.0;		//两个角点之间的最小欧式距离
    int block_size = 3;		//领域窗口大小
    bool Harris_use = false;		//是否适用harris
    double k = 0.06;		//经验系数（0.04-0.06）

    cv::Mat result_image = src.clone();
    cv::goodFeaturesToTrack(gray_src, corners, num_corners, qualityLevel, min_distance, cv::Mat(), block_size, Harris_use, k);

    for (int i = 0; i < corners.size(); i++)
        circle(result_image, corners[i], 15, (0, 0, 255), 3, 8, 0);

    cv::imshow(output_image, result_image);
}

int thre = 40;

void trackBar(int, void*);
void trackBar(int, void*)
{
//    std::vector<cv::KeyPoint> keypoints;
//    cv::Mat dst = gray_src.clone();
//    cv::Ptr<cv::FastFeatureDetector> detector = cv::FastFeatureDetector::create(thre);
//    detector->detect(gray_src,keypoints);
//    drawKeypoints(dst, keypoints, dst, cv::Scalar::all(-1), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
//    cv::imshow("output", dst);
}

void rectangle(cv::Mat img)
{
    if(img.empty())
        return;
    cv::Mat threasholdImg = img.clone();
    cv::threshold(threasholdImg,threasholdImg,40,255,cv::THRESH_BINARY);
    cv::Mat kenel = cv::getStructuringElement(cv::MORPH_RECT,cv::Size(5,5));
    cv::erode(threasholdImg,threasholdImg,kenel);
    //cv::GaussianBlur(img,img,cv::)
    //cv::Canny(img,img,100,200);
    cv::namedWindow("img2",cv::WINDOW_NORMAL);
    cv::imshow("img2",threasholdImg);
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarcy;

    cv::findContours(threasholdImg,contours,hierarcy,cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE);
    std::vector<cv::RotatedRect> box(contours.size());
    std::vector<cv::Rect> boundRect(contours.size());
    cv::Point2f rect[4];
    cv::Point2f rectMid[4];
    int xMin = threasholdImg.cols;
    int xMax = 0;
    int yMin = threasholdImg.rows;
    int yMax = 0;
    cv::Mat leftRoi;
    cv::Mat rightRoi;
    cv::Mat topRoi;
    cv::Mat bottomRoi;
    int radis = 80;
    std::vector<cv::Point> min_rectangle;
    std::vector<std::vector<cv::Point> > hull(contours.size());
    return;
    for(int i= 0;i < contours.size() ;i++)
    {
         box[i] = cv::minAreaRect(contours[i]);

         boundRect[i] = cv::boundingRect(contours[i]);
         if(boundRect[i].x <= 0 && boundRect[i].y <= 0&& boundRect[i].x >= img.rows && boundRect[i].y >= img.cols)
             continue;
         double area = cv::contourArea(contours[i]);

         if(area <= 1000)
             continue;
         qDebug() << "contours" << contours.size();
         float ratio =  (float)box[i].size.width /(float) box[i].size.height ;
         if( ratio < 0.85 && ratio > 1.15)
             continue;
         box[i].points(rect);
         qDebug() << "width" << "height" << box[i].size.width << box[i].size.height;
         qDebug() << "contoursArea" << area;
         //cv::rectangle(img,cv::Point(boundRect[i].x, boundRect[i].y), cv::Point(boundRect[i].x + boundRect[i].width, boundRect[i].y + boundRect[i].height), cv::Scalar(0, 255, 0), 2, 8);
         for(int j=0; j<4; j++)
        {
            //qDebug() << rect[j].x << rect[j].y;
            //rectMid[j] = (rect[j] + rect[(j+1)%4]) / 2;
            min_rectangle.push_back(rectMid[j]);
            rectMid[j].x = rectMid[j].x - radis;
            rectMid[j].y = rectMid[j].y - radis;
            cv::Mat src = img.clone();

            //cv::rectangle(img,cv::Rect(rectMid[j].x,rectMid[j].y,80,80),cv::Scalar(0,255,255),4);
            if(xMin > rectMid[j].x)
            {
                xMin = rectMid[j].x;
                leftRoi = src(cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2));
                cv::rectangle(img,cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2),cv::Scalar(255,255,255));
            }
            else if(xMax < rectMid[j].x)
            {
                xMax = rectMid[j].x;
                rightRoi = src(cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2));
                cv::rectangle(img,cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2),cv::Scalar(255,255,255));
            }else if(yMin > rectMid[j].y)
            {
                yMin = rectMid[j].y;
                topRoi = src(cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2));
                cv::rectangle(img,cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2),cv::Scalar(255,255,255));
            }else if(yMax < rectMid[j].y)
            {
                yMax = rectMid[j].y;
                bottomRoi = src(cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2));
                cv::rectangle(img,cv::Rect(rectMid[j].x,rectMid[j].y,radis * 2,radis * 2),cv::Scalar(255,255,255));
            }

            //SFRCalculation(roi,1);
            //cv::line(img, rect[j], rect[(j+1)%4], cv::Scalar(0, 0, 255), 8, 8);
         }
         double recArea = cv::contourArea(min_rectangle);
         qDebug() << recArea;
         //SFRCalculation(leftRoi,1);
         if(!rightRoi.empty())
         {
             cv::imshow("right",rightRoi);
             cv::Mat c = rightRoi.clone();
             cv::imshow("c",c);

              //qDebug() << "right" <<SFRCalculation(rightRoi,1);
             //SFRCalculation(rightRoi,1);
         }
         if(!leftRoi.empty())
         {
             cv::flip(leftRoi,leftRoi,-1);
             cv::imshow("left",leftRoi);
             //qDebug() << "left"<<SFRCalculation(leftRoi,1);
         }
         if(!topRoi.empty())
         {
             transpose(topRoi, topRoi);
             flip(topRoi, topRoi, 0);
             cv::imshow("top",topRoi);
             //qDebug() << "top" << SFRCalculation(topRoi,1);

         }
         if(!bottomRoi.empty())
         {
             cv::flip(bottomRoi,bottomRoi,0);
             cv::imshow("bottom",bottomRoi);
         }

    }


    cv::namedWindow("img",cv::WINDOW_NORMAL);
    cv::imshow("img",img);
}
void findPoint(cv::Point  &p, int flag, cv::Mat m, int *count)//点，flag，图，
{
    bool ok = true;
    int rows = m.rows / 2;
    int cols = m.cols / 2;
    switch (flag)
    {
    case 0:
    {
        while (ok && cols > 0)
        {
            if (m.at<uchar>(rows, cols) == 0)//255表示白点
            {
                p.x = rows;
                p.y = cols;
                //cout << "x=" << rows << ",y=" << cols << endl;
                ok = false;
            }
            cols--;
        }
        break;
    }
    case 1:
    {
        while (ok && cols < m.cols)
        {
            if (m.at<uchar>(rows, cols) == 0)
            {
                p.x = rows;
                p.y = cols;
                //cout << "x=" << rows << ",y=" << cols << endl;
                ok = false;
            }
            cols++;
        }
        break;
    }
    case 2:
    {
        while (ok && rows > 0 )
        {
            if (m.at<uchar>(rows, cols) == 0)
            {
                p.x = rows;
                p.y = cols;
                ok = false;
            }
            rows--;
        }
        break;
    }
    case 3:
    {
        while (ok && rows < m.rows)
        {
            if (m.at<uchar>(rows, cols) == 0)
            {
                p.x = rows;
                p.y = cols;

                //cout << "x=" << rows << ",y=" << cols << endl;
                ok = false;
            }
            rows++;
        }
        break;
    }
    default:
        break;
    }
}
void autolinkRect()
{
//    cv::Mat m = imageSource.clone();
//    int withd, height, Xoffset, Yoffset;
//    int rows = m.rows;
//    int cols = m.cols;
//    int R = sqrt(rows*rows + cols * cols) / 2;
//    double r_R = (rows / 2.000) / R;
//    double c_R = (cols / 2.000) / R;
//    cv::Point p = cv::Point(0, 0);
//    LinkCvPoint linkP;
//    withd = AreaW;
//    height = AreaH;

//    //0视场
//    p.x = cols / 2;
//    p.y = rows / 2;
//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(0, linkP);

//    //0.3视场
//    //左
//    p.x = cols / 2 - round(R*0.3);
//    p.y = rows / 2;
//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(1, linkP);
//    //右
//    p.x = cols / 2 + round(R*0.3);
//    p.y = rows / 2;
//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(2, linkP);
//    //上
//    p.x = cols / 2;
//    p.y = rows / 2 - round(R*0.3);
//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(3, linkP);
//    //下
//    p.x = cols / 2;
//    p.y = rows / 2 + round(R*0.3);
//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(4, linkP);

//    //0.7视场
//     //左上
//    p.x = cols / 2 - round(R*0.7*c_R);
//    p.y = rows / 2 - round(R*0.7*r_R);

//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;

//    setLinkCvRect(5, linkP);


//    p.x = cols / 2 + round(R*0.7*c_R);
//    p.y = rows / 2 - round(R*0.7*r_R);

//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(6, linkP);

//    //左下
//    p.x = cols / 2 - round(R*0.7*c_R);
//    p.y = rows / 2 + round(R*0.7*r_R);

//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(7, linkP);
//    p.x = cols / 2 + round(R*0.7*c_R);
//    p.y = rows / 2 + round(R*0.7*r_R);

//    Xoffset = p.y - height / 2;
//    Yoffset = p.x - withd / 2;
//    linkP.p = p;
//    linkP.withd = withd;
//    linkP.height = height;
//    linkP.Xoffset = Xoffset;
//    linkP.Yoffset = Yoffset;
//    setLinkCvRect(8, linkP);

}

cv::Mat src_img;
cv::Mat src_gray;
int thresh = 0;
int max_thresh = 255;
cv::RNG rng(12345);

void  RGBToGray(cv::Mat &src, cv::Mat &src_gray);
/**
 * @brief sfrCalculation
 * @param src 输入图片
 * @param hSfrVec 横向 SFR值
 * @param hEsfVec 横向 ESF值
 * @param hLsfVec 横向 LSF值
 * @param vSfrVec 纵向 SFR值
 * @param vEsfVec 纵向 ESF值

 * @param vLsfVec 纵向 LSF值
 * @param rate Gama值(越大对图像会越淡)
 * @param rectangleSize (横向与纵向正方矩形长宽大小,提高大小会提高SFR计算准度,一般要大于40像素，且是偶数）
 * @param showFlag 窗口显示
 */
void sfrCalculation(cv::Mat src,std::vector<double>&hSfrVec,std::vector<double>&hEsfVec,std::vector<double>&hLsfVec,std::vector<double>&vSfrVec,std::vector<double>&vEsfVec,std::vector<double>&vLsfVec,double rate,int rectangleSize,bool showFlag)
{
    if(src.empty())
        return;
    if(showFlag)
    {
        namedWindow("source", cv::WINDOW_AUTOSIZE);
        imshow("source", src);
    }

    cv::Mat threshold_output;
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::Mat src_gray;
    RGBToGray(src, src_gray);

    cv::Mat out_gray = src_gray.clone();

    /// 使用Threshold检测边缘
    cv::threshold(src_gray, threshold_output, 0, 255, cv::THRESH_BINARY_INV|cv::THRESH_OTSU);
    cv::Mat ker = cv::getStructuringElement(cv::MORPH_RECT,cv::Size(3,3));
    cv::erode(threshold_output,threshold_output,ker);
    if(showFlag)
    {
        cv::namedWindow("ThresholdOutput",cv::WINDOW_NORMAL);
        cv::imshow("ThresholdOutput",threshold_output);
    }


    cv::findContours(threshold_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0));
    std::vector<std::vector<cv::Point>> contours_poly(contours.size());
    std::vector<cv::Rect> boundRect(contours.size());
    bool flag = false;
    bool contoursSizeFlag = false;
    int count = 0;
    int currentPlotNum = 0;
    cv::Point p;
    cv::Point p2;
     cv::Mat temp = threshold_output.clone();
    cv::Mat drawing = temp.clone();
    hEsfVec.clear();
    hSfrVec.clear();
    hLsfVec.clear();
    vEsfVec.clear();
    vSfrVec.clear();
    vLsfVec.clear();
    if(contours.size() != 0)
    {
        if(contours.size() == 1)
            contoursSizeFlag = true;
        for (int i = 0; i < contours.size() ; i++)
        {
            qDebug()<<"contours.size()"<<contours.size();
            if(flag)
                continue;
            cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            cv::approxPolyDP(cv::Mat(contours[i]), contours_poly[i], 3, true);
            double minArea = cv::contourArea(contours_poly[i]);

            if(minArea <3000 || minArea >= 800000 )
                continue;

            boundRect[i] = boundingRect(cv::Mat(contours_poly[i]));

            if(showFlag)
                qDebug() << "boundRect"<<boundRect[i].x << temp.cols /3 << boundRect[i].y << temp.rows;
            if(boundRect[i].x > temp.cols / 3 && !contoursSizeFlag)
                continue;
            if((float)boundRect[i].width /(float) boundRect[i].height >2 || (float)boundRect[i].width / (float)boundRect[i].height <= 0.5)
                continue;
            if(showFlag)
            {
                qDebug() << "minArea" << minArea << (float)boundRect[i].width  / (float)boundRect[i].height;//645152
                cv::imshow(QString("temp%1").arg(currentPlotNum).toLocal8Bit().data(),temp(boundRect[i]));
            }

            findPoint(p,1,temp(boundRect[i]),&count);
            findPoint(p2,2,temp(boundRect[i]),&count);

            cv::Rect tempRect = cv::Rect(boundRect[i].x + p.y - rectangleSize/2,boundRect[i].y + p.x - rectangleSize/2,rectangleSize,rectangleSize);
            cv::Rect tempRect2 = cv::Rect(boundRect[i].x + p2.y - rectangleSize/2,boundRect[i].y + p2.x - rectangleSize/2,rectangleSize,rectangleSize);
            if(showFlag)
                qDebug() << "tempRect" << tempRect.x <<threshold_output.cols ;
            if( tempRect.y > 0 && tempRect.x > 0 && tempRect.x < threshold_output.cols  && tempRect.y < threshold_output.rows && tempRect.x + tempRect.width < threshold_output.cols && tempRect.y + tempRect.height < threshold_output.rows)
            {
                flag = true;
                cv::Mat tm = out_gray(tempRect);
                cv::rectangle(drawing, tempRect, color, 2, 8, 0);
                if(showFlag)
                    cv::imshow(QString("tm%1").arg(currentPlotNum).toLocal8Bit().data(),tm);

                SFRCalculation(tm, rate,hSfrVec,hEsfVec,hLsfVec);
                if(showFlag)
                {
                    QColor qc = QColor::fromHsl(rand()%360,rand()%256,rand()%200);
                    custplotesf->addGraph();
                    custplotesf->graph(currentPlotNum)->setPen(QPen(2));
                    custplotesf->graph(currentPlotNum)->setName(QString::number(hEsfVec[0]));
                    qDebug()<<"custplotesf->graph = "<<QString::number(hEsfVec[0]);
                    for(int i = 0;i < hEsfVec.size() ;i++)
                    {
                        if(custplotesf)
                        {
                            custplotesf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotesf->yAxis->rescale(true);

                            custplotesf->graph(currentPlotNum)->addData((double)i,hEsfVec[i]);

                        }
                    }

                    custplotesf->replot();

                    custplot->addGraph();
                    custplot->graph(currentPlotNum)->setPen(QPen(2));
                    custplot->graph(currentPlotNum)->setName(QString::number(hSfrVec[0]));

                    for(int i =0;i < hSfrVec.size();i++)
                    {

                        if(custplot)
                        {
                            custplot->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplot->yAxis->rescale(true);
                            custplot->graph(currentPlotNum)->addData((double)i/(double)hSfrVec.size(),hSfrVec[i]);

                        }
                    }
                    custplot->replot();

                    custplotlsf->addGraph();
                    custplotlsf->graph(currentPlotNum)->setPen(QPen(2));
                    custplotlsf->graph(currentPlotNum)->setName(QString::number(hLsfVec[0]));

                    for(int i = 0;i < hLsfVec.size() ;i++)
                    {
                        if(custplotlsf)
                        {
                            custplotlsf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotlsf->yAxis->rescale(true);
                            custplotlsf->graph(currentPlotNum)->addData((double)i,hLsfVec[i]);
                        }
                    }
                    custplotlsf->replot();
                    currentPlotNum++;
                }



            }

            if( tempRect2.y > 0 && tempRect2.x > 0 && tempRect2.x < threshold_output.cols && tempRect2.y < threshold_output.rows && tempRect2.x + tempRect2.width < threshold_output.cols && tempRect2.y + tempRect2.height < threshold_output.rows)
            {
                cv::Mat tm = out_gray(tempRect2);
                cv::rectangle(drawing, tempRect2, color, 2, 8, 0);
                transpose(tm, tm);// 翻转模式，flipCode == 0垂直翻转（沿X轴翻转），flipCode>0水平翻转（沿Y轴翻转），flipCode<0水平垂直翻转（先沿X轴翻转，再沿Y轴翻转，等价于旋转180°）
                flip(tm, tm, 1);
                SFRCalculation(tm, rate,vSfrVec,vEsfVec,vLsfVec);
                if(showFlag)
                {
                    custplot->addGraph();
                    custplot->graph(currentPlotNum)->setPen(QPen(1));
                    //custplot->graph(currentPlotNum)->setName(QString::number(vSfrVec[0]));
                    custplot->graph(currentPlotNum)->setName(QString::number(66666));
                    for(int i =0;i < vSfrVec.size();i++)
                    {

                        if(custplot)
                        {
                            custplot->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplot->yAxis->rescale(true);
                            custplot->graph(currentPlotNum)->addData((double)i/(double)vSfrVec.size(),vSfrVec[i]);

                        }
                    }
                    custplot->replot();

                    custplotesf->addGraph();
                    custplotesf->graph(currentPlotNum)->setPen(QPen(1));
                    //custplotesf->graph(currentPlotNum)->setName(QString::number(vEsfVec[0]));
                    custplotesf->graph(currentPlotNum)->setName(QString::number(77777));
                    for(int i = 0;i < vEsfVec.size() ;i++)
                    {
                        if(custplotesf)
                        {
                            custplotesf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotesf->yAxis->rescale(true);

                            custplotesf->graph(currentPlotNum)->addData((double)i,vEsfVec[i]);

                        }
                    }
                    custplotesf->replot();

                    custplotlsf->addGraph();
                    custplotlsf->graph(currentPlotNum)->setPen(QPen(1));
                    //custplotlsf->graph(currentPlotNum)->setName(QString::number(vLsfVec[0]));
                    custplotlsf->graph(currentPlotNum)->setName(QString::number(88888));

                    for(int i = 0;i < vLsfVec.size() ;i++)
                    {
                        if(custplotlsf)
                        {
                            custplotlsf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotlsf->yAxis->rescale(true);
                            custplotlsf->graph(currentPlotNum)->addData((double)i,vLsfVec[i]);
                        }
                    }
                    custplotlsf->replot();
                    currentPlotNum++;
                }

                if(showFlag)
                    cv::imshow("tm",tm);
            }

            //minEnclosingCircle(contours_poly[i], center[i], radius[i]);
        }

    }
    if(showFlag)
    {
        cv::imshow("drawing",drawing);
        cv::waitKey(0);
    }

}
//#define ok
#ifdef ok
/// 函数声明
void thresh_callback(int, void*);


/** @thresh_callback 函数 */
void thresh_callback(int val, void*)
{
        return;

        cv::Mat threshold_output;
        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::Mat out_gray = src_gray.clone();
        /// 使用Threshold检测边缘
        cv::threshold(src_gray, threshold_output, thresh, 255, cv::THRESH_BINARY_INV|cv::THRESH_OTSU);
        cv::Mat threshold_temp = threshold_output.clone();
        cv::Mat ker = cv::getStructuringElement(cv::MORPH_RECT,cv::Size(3,3));
        cv::erode(threshold_output,threshold_output,ker);
        /// 显示到窗口中
        namedWindow("threshold_output", cv::WINDOW_NORMAL);
        imshow("threshold_output", threshold_output);

        /// 找到轮廓
        cv::findContours(threshold_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0));

        /// 计算矩
        std::vector<cv::Moments> mu(contours.size());
        for (int i = 0; i < contours.size(); i++)
        {
            mu[i] = moments(contours[i], false);
        }

        ///  计算中心矩:
        std::vector<cv::Point2f> mc(contours.size());
        for (int i = 0; i < contours.size(); i++)
        {
            mc[i] = cv::Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
        }
        cv::Point p;
        cv::Point p2;
        int count = 0;
        /// 绘制轮廓
        cv::Mat drawing = cv::Mat::zeros(src_gray.size(), CV_8UC3);
        for (int i = 0; i< contours.size(); i++)
        {
            cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, cv::Point());
            circle(drawing, mc[i], 4, color, -1, 8, 0);

        }

        /// 显示到窗口中
        cv::namedWindow("Contours", cv::WINDOW_NORMAL);
        cv::imshow("Contours", drawing);

        /// 通过m00计算轮廓面积并且和OpenCV函数比较
        printf("\t Info: Area and Contour Length \n");
        for (int i = 0; i< contours.size(); i++)
        {
            //printf(" * Contour[%d] - Area (M_00) = %.2f - Area OpenCV: %.2f - Length: %.2f \n", i, mu[i].m00, contourArea(contours[i]), arcLength(contours[i], true));
            cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, cv::Point());
            circle(drawing, mc[i], 4, color, -1, 8, 0);
        }
        /// 多边形逼近轮廓 + 获取矩形和圆形边界框
        std::vector<std::vector<cv::Point>> contours_poly(contours.size());
        std::vector<cv::Rect> boundRect(contours.size());
        std::vector<cv::Point2f>center(contours.size());
        std::vector<float>radius(contours.size());
        cv::Mat temp = threshold_temp.clone();
        std::vector<double> res;
        std::vector<double> esfVec;
        std::vector<double> lsfVec;
        int currentPlotNum = 0;
        qDebug() << ">>>>>>>>>"<<contours.size();
        bool flag = false;
        if(contours.size() != 0)
        {
            for (int i = 0; i < contours.size() ; i++)
            {
                if(flag)
                    continue;
                cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
                approxPolyDP(cv::Mat(contours[i]), contours_poly[i], 3, true);
                double minArea = cv::contourArea(contours_poly[i]);

                if(minArea <3000 || minArea >= 800000 )
                    continue;
                res.clear();
                esfVec.clear();
                lsfVec.clear();
                boundRect[i] = boundingRect(cv::Mat(contours_poly[i]));

                //if(showFlag)
                    qDebug() << "boundRect"<<boundRect[i].x << temp.cols /3 << boundRect[i].y << temp.rows;
                if(boundRect[i].x > temp.cols / 3)
                    continue;
                if((float)boundRect[i].width /(float) boundRect[i].height >2 || (float)boundRect[i].width / (float)boundRect[i].height <= 0.5)
                    continue;
                //if(showFlag)
                    qDebug() << "minArea" << minArea << (float)boundRect[i].width  / (float)boundRect[i].height;//645152
                cv::imshow(QString("temp%1").arg(currentPlotNum).toLocal8Bit().data(),temp(boundRect[i]));

                findPoint(p,1,temp(boundRect[i]),&count);
                findPoint(p2,2,temp(boundRect[i]),&count);

                cv::Rect tempRect = cv::Rect(boundRect[i].x + p.y - 20,boundRect[i].y + p.x - 20,40,40);
                cv::Rect tempRect2 = cv::Rect(boundRect[i].x + p2.y - 20,boundRect[i].y + p2.x - 20,40,40);
                //if(showFlag)
                    qDebug() << "tempRect" << tempRect.x <<threshold_output.cols ;
                if( tempRect.y > 0 && tempRect.x > 0 && tempRect.x < threshold_output.cols  && tempRect.y < threshold_output.rows && tempRect.x + tempRect.width < threshold_output.cols && tempRect.y + tempRect.height < threshold_output.rows)
                {
                    flag = true;
                    cv::Mat tm = out_gray(tempRect);
                    cv::rectangle(drawing, tempRect, color, 2, 8, 0);
                   // if(showFlag)
                        cv::imshow(QString("tm%1").arg(currentPlotNum).toLocal8Bit().data(),tm);

                    SFRCalculation(tm, 1,res,esfVec,lsfVec);

                    QColor qc = QColor::fromHsl(rand()%360,rand()%256,rand()%200);
                    custplotesf->addGraph();
                    custplotesf->graph(currentPlotNum)->setPen(QPen(qc));
                    custplotesf->graph(currentPlotNum)->setName(QString::number(esfVec[0]));
                    for(int i = 0;i < esfVec.size() ;i++)
                    {
                        if(custplotesf)
                        {
                            custplotesf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotesf->yAxis->rescale(true);

                            custplotesf->graph(currentPlotNum)->addData((double)i,esfVec[i]);

                        }
                    }

                    custplotesf->replot();

                    custplot->addGraph();
                    custplot->graph(currentPlotNum)->setPen(QPen(qc));
                    custplot->graph(currentPlotNum)->setName(QString::number(res[0]));

                    for(int i =0;i < res.size();i++)
                    {

                        if(custplot)
                        {
                            custplot->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplot->yAxis->rescale(true);
                            custplot->graph(currentPlotNum)->addData((double)i/(double)res.size(),res[i]);

                        }
                    }
                    custplotlsf->addGraph();
                    custplotlsf->graph(currentPlotNum)->setPen(QPen(qc));
                    custplotlsf->graph(currentPlotNum)->setName(QString::number(lsfVec[0]));

                    for(int i = 0;i < lsfVec.size() ;i++)
                    {
                        if(custplotlsf)
                        {
                            custplotlsf->xAxis->rescale(true);//调整X轴的范围，使之能显示出所有的曲线的X值
                            custplotlsf->yAxis->rescale(true);
                            custplotlsf->graph(0)->addData((double)i,lsfVec[i]);
                        }
                    }
                    custplotlsf->replot();
                    currentPlotNum++;

                }
                if( tempRect2.y > 0 && tempRect2.x > 0 && tempRect2.x < threshold_output.cols && tempRect2.y < threshold_output.rows && tempRect2.x + tempRect2.width < threshold_output.cols && tempRect2.y + tempRect2.height < threshold_output.rows)
                {
                    cv::Mat tm = out_gray(tempRect2);
                    cv::rectangle(drawing, tempRect2, color, 2, 8, 0);
                    transpose(tm, tm);// 翻转模式，flipCode == 0垂直翻转（沿X轴翻转），flipCode>0水平翻转（沿Y轴翻转），flipCode<0水平垂直翻转（先沿X轴翻转，再沿Y轴翻转，等价于旋转180°）
                    flip(tm, tm, 1);
                    //cv::imshow("tm2" + std::to_string(i),tm);
                    //std::cout <<">>>>>>>>>>>>y  "<< i <<SFRCalculation(tm, 1) << std::endl;
                }

                minEnclosingCircle(contours_poly[i], center[i], radius[i]);
            }

        }

        /// 画多边形轮廓 + 包围的矩形框 + 圆形框

        for (int i = 0; i< contours.size(); i++)
        {
            cv::Scalar color = cv::Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            //drawContours(drawing, contours_poly, i, color, 1, 8, std::vector<cv::Vec4i>(), 0, cv::Point());
            cv::rectangle(drawing, boundRect[i].tl(), boundRect[i].br(), color, 2, 8, 0);
            //circle(drawing, center[i], (int)radius[i], color, 2, 8, 0);
        }

        /// 显示在一个窗口
        cv::namedWindow("Contours", cv::WINDOW_NORMAL);
        cv::imshow("Contours", drawing);


}
#endif

void initCustomPlot()
{
    custplot = new QCustomPlot();
    custplot->show();
    custplot->resize(300,300);
    custplot->xAxis->setAutoTickStep(true);
    custplot->xAxis->setTickStep(5);
    custplot->yAxis->setAutoTickStep(true);
    custplot->axisRect()->setupFullAxesBox();
    custplot->xAxis->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);//设置图例位置左上角
    custplot->legend->setVisible(true);//图例显示与否
    custplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);


    custplotesf = new QCustomPlot();
    custplotesf->show();
    custplotesf->resize(300,300);

    custplotesf->xAxis->setAutoTickStep(true);
    custplotesf->xAxis->setTickStep(5);
    custplotesf->yAxis->setAutoTickStep(true);
    custplotesf->axisRect()->setupFullAxesBox();
    custplotesf->xAxis->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);//设置图例位置左上角
    custplotesf->legend->setVisible(true);//图例显示与否
    custplotesf->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);


    custplotlsf = new QCustomPlot();
    custplotlsf->show();
    custplotlsf->resize(300,300);
    custplotlsf->addGraph();
    custplotlsf->graph(0)->setPen(QPen(Qt::red));
    custplotlsf->graph(0)->setName("MTF");
    custplotlsf->xAxis->setAutoTickStep(true);
    custplotlsf->xAxis->setTickStep(5);
    custplotlsf->yAxis->setAutoTickStep(true);
    custplotlsf->axisRect()->setupFullAxesBox();
    custplotlsf->xAxis->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);//设置图例位置左上角
    custplotlsf->legend->setVisible(true);//图例显示与否
    custplotlsf->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

}


void RGBToGray(cv::Mat &src, cv::Mat &des)
{
    //Gray = (R*30 + G*59 +B*11 +50)/100;
    des.create(src.size(), CV_8UC1);
    for (int r = 0; r < src.rows; r++)
    {
        for (int c = 0; c < src.cols; c++)
        {
            cv::Vec3b &m = src.at<cv::Vec3b>(r, c);
            int gray = (m[2] * 0 + m[1] * 100 + m[0] * 0 + 50) / 100;
            des.at<uchar>(r, c) = gray;
        }
    }

}
/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[])
{
        QApplication a(argc,argv);
        qDebug() << "Main";
        initCustomPlot();
    /// 载入原图像, 返回3通道图像
        src_img = cv::imread("D:/QTcode/SFRTest/SFR/SFRimage/r3.jpg", 1);
        if(src_img.empty())
            qDebug() << "Image is Empty!";
        std::vector<double> sfrVec;
        std::vector<double> esfVec;
        std::vector<double> lsfVec;
        std::vector<double> vsfrVec;
        std::vector<double> vesfVec;
        std::vector<double> vlsfVec;
        sfrCalculation(src_img,sfrVec,esfVec,lsfVec,vsfrVec,vesfVec,vlsfVec,1,40,true);
        cv::waitKey(0);
        return -1;

}

