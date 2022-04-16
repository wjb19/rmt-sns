// WJB 2020,2021
//
// take output grid and map satellite image pixels to web mercator projection
//
//
// requirements eg., on ubuntu:
// sudo apt install libboost-all-dev
// sudo apt install libeigen3-dev
// also needs armadillo, knn-cpp and pugixml
//


// openstreet map tile range no zoom
// (-180, +85.0511) - (+180, -85.0511)

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <iostream>
#include <unordered_map>
#include <armadillo>


#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std::chrono;


/// https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames#Lon..2Flat._to_bbox
int long2tilex(double lon, int z) 
{ 
  return (int)((lon + 180.0) / 360.0 * (1<<z)); 
}

int lat2tiley(double lat, int z)
{ 
  double latrad = lat * M_PI/180.0;
  return (int)((1.0 - log(tan(latrad) + 1/cos(latrad)) / M_PI) / 2.0 * (1<<z)); 
}

double tilex2long(double x, int z) 
{
  return x / (double)(1 << z) * 360.0 - 180;
}

double tiley2lat(double y, int z) 
{
  double n = M_PI - 2.0 * M_PI * y / (double)(1 << z);
  return 180.0 / M_PI * atan(0.5 * (exp(n) - exp(-n)));
}



struct lonLat
{
  double lat,lon;
  double x,y,z;
  double imx, imy;
  bool inview;

  float pvalue;

  int tilex,tiley;
  int tilexp,tileyp;


  lonLat(std::vector<double>& in, int zoom)
  {
    //ret(index,0) = 2.;
    //ret(index,1) = (double) subtype;
    //ret(index,2) = p.lon;
    //ret(index,3) = p.lat;
    //ret(index,4) = p.x;
    //ret(index,5) = p.y;
    //ret(index,6) = p.z;
    //ret(index,7) = p.imx;
    //ret(index,8) = p.imy;

    lon =in[2];
    lat =in[3];
    x   =in[4];
    y   =in[5];
    z   =in[6];
    imx =in[7];
    imy =in[8];

    tilex = long2tilex(lon,zoom);
    tiley = lat2tiley(lat, zoom);

    zoom+=8;

    tilexp = long2tilex(lon,zoom);
    tileyp = lat2tiley(lat, zoom);

  }

  bool operator==(const lonLat& ll) const
  {
    return ll.lat == lat && ll.lon == lon;
  }
};


struct xyHash
{
  size_t operator()(const std::pair<int,int>& xy) const
  { 
    return (std::hash<int>()(xy.first)) ^  
      (std::hash<int>()(xy.second)); 
  } 
};

struct llHash
{
  size_t operator()(const std::pair<double,double>& ll) const
  { 
    return (std::hash<double>()(ll.first)) ^  
      (std::hash<double>()(ll.second)); 
  } 
};

///
arma::Mat<float> createTile(int x, int y, std::vector<lonLat>& pixels, arma::Mat<float>& image)
{

  arma::Mat<float> ret(256,256,arma::fill::ones),w(256,256,arma::fill::ones);

  double ix = image.n_rows;
  double iy = image.n_cols;

  ret*=255.;

  for (auto &p : pixels)
  {
    int i = p.tilexp - p.tilex*256;
    int j = p.tileyp - p.tiley*256;

    int li =  static_cast<int>(ix + p.imx - ix/2.);  
    int lj =  static_cast<int>(iy + p.imy - iy/2.);  

    ret(j,i)+=image(li,lj);

    w(j,i)++;
  }

  for (int i=0; i<256; i++)
    for (int j=0; j<256; j++)
      ret(j,i)/=w(j,i);

  return ret;

}



using namespace boost::program_options;

int main(int argc, char * argv[])
{

  int z=1;

  options_description desc("webMerc - create web mercator projection from satellite image");
  desc.add_options()
    ("help", "produce help message")
    ("z", value<int>(), "zoom")
    ("in", value<std::string>(), "input image filename (required)")
    ("out", value<std::string>(), "output tiles basename (default=./)");

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  std::string base="./";

  if (vm.count("help"))
  {
    std::cerr << desc << std::endl;
    return 0;
  }

  if (vm.count("out"))
    base = vm["out"].as<std::string>();

  if (vm.count("z"))
    z = vm["z"].as<int>();

  std::string input_image;
  if (vm.count("in"))
    input_image = vm["in"].as<std::string>();
  else
  {
    std::cerr << desc << std::endl;
    return 0;
  }

  std::vector<double> buf(9);
  std::vector<lonLat> grid;

  std::unordered_map<std::pair<int,int> , std::vector<lonLat> ,xyHash> tile_pixels;
  std::unordered_map<std::pair<int,int> , std::vector<lonLat> ,xyHash> tile_pixels_complete;
  std::unordered_map<std::pair<double,double>, lonLat, llHash> ll_pixel;


  size_t num_points=0; 

  double mnx=1e20,mxx=-1e20,mny=1e20,mxy=-1-20;
  while(std::cin)
  {
    std::cin.read((char*)buf.data(),sizeof(double)*9);
    lonLat ll(buf,z);
    std::pair<int,int> key1(ll.tilex,ll.tiley);
    std::pair<double,double> key2(ll.lon,ll.lat);
    tile_pixels[key1].push_back(ll);
    ll_pixel.insert(std::pair<std::pair<double,double> , lonLat>(key2,ll));
    if (ll.imx > mxx) mxx=ll.imx;
    if (ll.imx < mnx) mnx=ll.imx;
    if (ll.imy > mxy) mxy=ll.imy;
    if (ll.imy < mny) mny=ll.imy;
    num_points++;
  }

  size_t index=0;
  for (auto it=ll_pixel.begin(); it!=ll_pixel.end(); it++)
  {
    lonLat ll = it->second;
    index++;
  }


  size_t thres = 32000;


  cv::Mat im = cv::imread(input_image.c_str(),cv::IMREAD_GRAYSCALE);

  arma::Mat<float> image(im.rows,im.cols);

  for (int i=0; i<im.rows; i++)
    for (int j=0; j<im.cols; j++)
      image(i,j) = static_cast<float>(im.at<unsigned char>(i,j));



  for (auto &n : tile_pixels)
  {
    if (n.second.size() > thres)
    {
      std::string filename = base+"/"+std::to_string(z)+"_";
      filename+=std::to_string(n.first.first)+"_";
      filename+=std::to_string(n.first.second)+".png";
      auto tile = createTile(n.first.first,n.first.second,n.second,image);

      cv::Mat out(256,256,CV_8UC1);

      for (int i=0; i<256; i++)
        for (int j=0; j<256; j++)
          out.at<unsigned char>(i,j) = static_cast<unsigned char>(tile(i,j));

      cv::imwrite(filename.c_str(),out);

    }
  } 

  return 0;
}

