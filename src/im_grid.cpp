// WJB 2020,2021
//
// work out what grid the satellite camera is looking at, so we can ultimately map pixels 
// correctly and create web mercator projection, suited to openstreet map
//
//
// references
// https://www.edn.com/digital-camera-design-determine-pixel-size-focal-ratio-and-field-of-view/
// https://courses.cs.washington.edu/courses/cse455/09wi/Lects/lect5.pdf
// https://en.wikipedia.org/wiki/Keyhole_Markup_Language
// https://www.ngs.noaa.gov/GRD/GPS/DOC/toc.html
//
// requirements eg., on ubuntu:
// sudo apt install libboost-all-dev
// sudo apt install libeigen3-dev
// also needs armadillo, knn-cpp and pugixml
//


#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <iostream>
#include <unordered_map>

#include "knn/kdtree_minkowski.h"
#include "pugixml.hpp"
#include "armadillo"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
typedef Eigen::MatrixXd Matrix;
typedef knn::Matrixi Matrixi;

#define aa 6378137.0
#define aa2  aa*aa
#define aa4 aa2*aa2
#define bb 6356752.3142
#define bb2 bb*bb
#define bb4 bb2*bb2

// create world x,y,z coordinates from lon/lat

inline void earthXYZ(double lon, double lat, 
    double& x, double& y, double& z)
{
  //lat = 90-lat;
  lat /= 180.;
  lat *= M_PI;

  //lon += 180;
  lon /= 180.;
  lon *= M_PI;

  double c = cos(lat);
  double s = sin(lat);
  double c2 = c*c;
  double s2 = s*s;

  double rad = sqrtf((aa4 * c2 + bb4 * s2) / (aa2 * c2 + bb2 * s2));

  s = sin(M_PI/2.-lat);
  c = cos(M_PI/2.-lat);

  x= rad* s * cos(lon+M_PI);
  y= rad* s * sin(lon+M_PI);

  z=rad* c;
};



struct lonLat
{
  double lat,lon;
  double x,y,z;
  double imx, imy;
  bool inview;

  lonLat(double ln, double lt): 
    lon(ln),lat(lt)
  {
    earthXYZ(lon,lat,x,y,z);
    inview=false;
  }

  bool operator==(const lonLat& ll) const
  {
    return ll.lat == lat && ll.lon == lon;
  }
};


struct lonLatHash
{
  size_t operator()(const lonLat& ll) const
  { 
    return (std::hash<double>()(ll.lat)) ^  
      (std::hash<double>()(ll.lon)); 
  } 
};



// create world frame to spacecraft frame rotation matrix

inline arma::Mat<double> w2c_rot(double lon, double lat, 
    double ya=0., double pa=0., double ra=0.)
{

  ya /= 180.;
  ya *= M_PI;

  pa /= 180.;
  pa *= 180.;

  ra /= 180.;
  ra *= M_PI;

  lonLat ll(lon,lat);

  // tait - bryan z-y'-x'' applied as sequence yaw - pitch - roll

  // FIRST step : use craft lon (initial yaw) /lat (initial pitch) position to orient new frame
  // check that :
  // *new z axis is parallel with ray through lon lat in world coords
  // *new x/y plane orthogonal to same ray
  //
  // SECOND step : apply craft yaw / pitch / roll relative to this frame

  lon += 180.;
  lon /= 180.;
  lon *= M_PI;

  lat = 90.-lat;
  lat /= 180.;
  lat *= M_PI;

  arma::Mat<double> yaw(3,3);

  yaw(0,0) = cos(lon);
  yaw(0,1) = sin(lon);
  yaw(0,2) = 0.;

  yaw(1,0) = -sin(lon);
  yaw(1,1) = cos(lon);
  yaw(1,2) = 0.;

  yaw(2,0) = 0.;
  yaw(2,1) = 0.;
  yaw(2,2) = 1.;

  arma::Mat<double> pit(3,3);

  pit(0,0) = cos(lat);
  pit(0,1) = 0.;
  pit(0,2) = -sin(lat);

  pit(1,0) = 0.;
  pit(1,1) = 1.;
  pit(1,2) = 0.;

  pit(2,0) = sin(lat);
  pit(2,1) = 0.;
  pit(2,2) = cos(lat);

  arma::Mat<double> rot = pit*yaw;
  arma::Mat<double> irot = inv(rot);
  arma::Mat<double> test_vector(3,1);

  // z parallel to ray in world

  test_vector(0,0) = 0.;
  test_vector(1,0) = 0.;
  test_vector(2,0) = 6400.*1000.;

  arma::Mat<double> o = irot*test_vector;

  double mag1 = sqrt(ll.x*ll.x + ll.y*ll.y + ll.z*ll.z);
  double mag2 = sqrt(o(0)*o(0) + o(1)*o(1) + o(2)*o(2));

  double ct = (ll.x*o(0)+ll.y*o(1)+ll.z*o(2)) / (mag1*mag2);

  assert(fabs(ct-1.)<1e-10);

  // x orthog to ray in world

  test_vector(0,0) = 0.;
  test_vector(1,0) = 6400.*1000.;
  test_vector(2,0) = 0.;

  o = irot*test_vector;
  mag2 = sqrt(o(0)*o(0) + o(1)*o(1) + o(2)*o(2));

  ct = (ll.x*o(0)+ll.y*o(1)+ll.z*o(2)) / (mag1*mag2);

  assert(fabs(ct)<1e-10);

  // y orthog to ray in world

  test_vector(0,0) = 6400.*1000.;
  test_vector(1,0) = 0.;
  test_vector(2,0) = 0.;

  o = irot*test_vector;
  mag2 = sqrt(o(0)*o(0) + o(1)*o(1) + o(2)*o(2));
  ct = (ll.x*o(0)+ll.y*o(1)+ll.z*o(2)) / (mag1*mag2);
  assert(fabs(ct)<1e-10);

  yaw(0,0) = cos(ya);
  yaw(0,1) = sin(ya);
  yaw(0,2) = 0.;

  yaw(1,0) = -sin(ya);
  yaw(1,1) = cos(ya);
  yaw(1,2) = 0.;

  yaw(2,0) = 0.;
  yaw(2,1) = 0.;
  yaw(2,2) = 1.;

  pit(0,0) = cos(pa);
  pit(0,1) = 0.;
  pit(0,2) = -sin(pa);

  pit(1,0) = 0.;
  pit(1,1) = 1.;
  pit(1,2) = 0.;

  pit(2,0) = sin(pa);
  pit(2,1) = 0.;
  pit(2,2) = cos(pa);

  arma::Mat<double> roll(3,3);

  roll(0,0) = 1.;
  roll(0,1) = 0.;
  roll(0,2) = 0.;

  roll(1,0) = 0.;
  roll(1,1) = cos(ra);
  roll(1,2) = sin(ra);

  roll(2,0) = 0.;
  roll(2,1) = -sin(ra);
  roll(2,2) = cos(ra);

  arma::Mat<double> output = roll*pit*yaw*rot;
  return output;
}


arma::Mat<double> inline camera_matrix(double lon, double lat, 
    double alt, double size, double yaw, double pitch, double roll)
{

  // euler rotation matrix
  arma::Mat<double> cmat = w2c_rot(lon, lat, yaw, pitch, roll);

  // intrinsic contribution
  arma::Mat<double> Mint(3,3); 

  Mint.zeros();
  Mint(0,0) = alt / (size);
  Mint(1,1) = alt / (size);
  Mint(2,2) = 1.;

  // extrinsic
  arma::Mat<double> Mext(3,4);
  // translation vector
  lat /= 180.;
  lat *= M_PI;

  lon /= 180.;
  lon *= M_PI;

  double c = cos(lat);
  double s = sin(lat);
  double c2 = c*c;
  double s2 = s*s;

  double rad = sqrtf((aa4 * c2 + bb4 * s2) / (aa2 * c2 + bb2 * s2)) + alt*1000.;

  s = sin(M_PI/2.-lat);
  c = cos(M_PI/2.-lat);

  double x= rad* s * cos(lon+M_PI);
  double y= rad* s * sin(lon+M_PI);
  double z= rad* c;

  Mext.submat(0,0,2,2) = cmat;

  Mext(0,3) = -cmat(0,0)*x - cmat(0,1)*y - cmat(0,2)*z;
  Mext(1,3) = -cmat(1,0)*x - cmat(1,1)*y - cmat(1,2)*z;
  Mext(2,3) = -cmat(2,0)*x - cmat(2,1)*y - cmat(2,2)*z;

  arma::Mat<double> M = Mint*Mext;

  return M;

}


typedef std::vector<lonLat> regions;
typedef std::vector<regions > all_regions;
typedef std::unordered_map<lonLat,std::string,lonLatHash> coords_country;
typedef std::unordered_map<std::string,all_regions > country_coords;
typedef std::unordered_map<std::string,regions > lon_lat_coords;

// create lookup tables for mapping between coordinates <-> country name
// from kml file of country polygons

inline int create_luts(std::string& input_world_polygons_file, 
    coords_country& coords_country_map, country_coords& country_coords_map)
{


  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(input_world_polygons_file.c_str());

  if (!result)
    return 0;

  std::cerr << "loaded world polys; parsing file..." << std::endl;

  size_t total_points = 0;
  auto folders = doc.child("kml").child("Document");
  std::string country;

  for (pugi::xml_node f = folders.child("Folder"); f; f = f.next_sibling("Folder"))
  {

    for (pugi::xml_node_iterator it= f.begin(); it!=f.end(); it++)
    {
      country = it->child_value("name");

      std::vector<std::string> all_polys;
      size_t num_all_polygons = 0;

      // single polys
      for (pugi::xml_node pl = it->child("Polygon"); pl; pl = pl.next_sibling("Polygon"))
      {
        auto c = pl.child("outerBoundaryIs").child("LinearRing");
        std::string coords =  c.child_value("coordinates");
        all_polys.push_back(coords);
        num_all_polygons++;
      }

      // multi polys
      auto mpl  = it->child("MultiGeometry");
      for (pugi::xml_node pl = mpl.child("Polygon"); pl; pl = pl.next_sibling("Polygon"))
      {
        auto c = pl.child("outerBoundaryIs").child("LinearRing");
        std::string coords =  c.child_value("coordinates");
        all_polys.push_back(coords);
        num_all_polygons++;
      }

      all_regions all;
      for (auto& x: all_polys)
      {

        std::istringstream split(x);
        std::vector<std::string> tokens;
        for (std::string each; std::getline(split, each, ' '); tokens.push_back(each));

        std::vector<lonLat> region;
        for (auto&cc : tokens)
        {

          std::istringstream split2(cc);
          std::vector<std::string> lonlatalt;
          for (std::string each; std::getline(split2, each, ','); lonlatalt.push_back(each));

          if (lonlatalt.size() != 3) 
            continue;

          total_points++;
          double lon = stod(lonlatalt[0]);
          double lat = stod(lonlatalt[1]);
          double x,y,z;

          lonLat ll(lon,lat);
          coords_country_map.insert(std::pair<lonLat,std::string>{ll,country});       

          //std::cerr << lon << " " << lat << " " << total_points << " " << coords_country_map.size() 
          //  << " " << total_points - coords_country_map.size() << std::endl;

          region.push_back(ll);
        }
        all.push_back(region);
      }

      //std::cerr << country << " " << all.size() << " " << num_all_polygons << std::endl;
      country_coords_map.insert(std::pair<std::string,all_regions>{country,all});
    }
  }

  std::cerr << "...done" << std::endl;
  return 1;
}

// generate a lat/lon grid based on camera eye location, maximum distance and resolution 

void gen_grid(double eye_lat, double eye_lon, 
    double distance, double resolution, lon_lat_coords& lons, lon_lat_coords& lats)
{

  double rad = 6400.*1000.;
  double ang = atan (distance / rad);

  double ang_deg  = (double) (((int)( ang / M_PI * 180.) / 10u +1u)* 10u);
  size_t points = (size_t) ceil(ang_deg / resolution);

  std::cerr << ang_deg << " " << resolution << " " << points << std::endl;

  double start_lat = round(eye_lat - ang_deg / 2.);
  double start_lon = round(eye_lon - ang_deg / 2.);

  double lon = start_lon-resolution;

  for (int i=0; i<points; i++)
  {
    lon += resolution;

    if (lon < -180) lon += 360;
    else if (lon >  180) lon -= 360;
    else{}

    std::vector<lonLat> ll;
    double lat = start_lat-resolution;
    std::string label = "lon_"+std::to_string(lon) + "_";
    for (int j=0; j<points; j++)
    {
      lat  += resolution;

      double sign = (lat < 0) ? -1. : 1.;
      lat = fabs(lat);
      if (lat > 90) lat = 180 - lat;

      lat *= sign;

      ll.emplace_back(lon,lat);

      if (j==0)
        label+=std::to_string(lat) + "_";
    }

    label+=std::to_string(lat);
    lons.insert({label,ll});
  }

  double lat = start_lat-resolution;

  for (int i=0; i<points; i++)
  {

    lat  += resolution;

    double sign = (lat < 0) ? -1. : 1.;
    lat = fabs(lat);
    if (lat > 90) lat = 180 - lat;

    lat *= sign;

    std::vector<lonLat> ll;
    double lon = start_lon-resolution;
    std::string label = "lat_"+std::to_string(lat) + "_";
    for (int j=0; j<points; j++)
    {
      lon += resolution;

      if (lon < -180) lon += 360;
      else if (lon >  180) lon -= 360;
      else{}

      ll.emplace_back(lon,lat);
      if (j==0)
        label+=std::to_string(lon) + "_";
    }

    label+=std::to_string(lon);
    lats.insert({label,ll});
  }
}

// given a location, find nearest neighbor countries

std::vector<std::string> inline find_nns(coords_country& map, 
    double eye_lat, double eye_lon, double distance)
{

  std::vector<std::string> ret;
  std::cerr << "creating tree for finding nn countries.." << std::endl;

  size_t num_points = map.size(), index=0;
  Matrix dataPoints(3, num_points);

  double x,y,z;
  earthXYZ(eye_lon,eye_lat,x,y,z);

  for (auto it=map.begin(); it!=map.end(); it++)
  {
    lonLat ll = it->first;

    dataPoints(0,index)=ll.x;
    dataPoints(1,index)=ll.y;
    dataPoints(2,index)=ll.z;

    index++;
  }

  knn::KDTreeMinkowski<double, knn::EuclideanDistance<double>> kdtree(dataPoints);

  kdtree.setBucketSize(16);
  kdtree.setCompact(false);
  kdtree.setBalanced(false);
  kdtree.setSorted(true);
  kdtree.setTakeRoot(true);
  kdtree.setMaxDistance(distance);
  kdtree.setThreads(2);
  kdtree.build();

  std::cerr << "..done" << std::endl;

  Matrix queryPoints(3, 1);
  queryPoints << x, y, z;

  Matrixi indices;
  Matrix distances;
  kdtree.query(queryPoints, 400, indices, distances);


  //std::cerr
    //<< "Data points:" << std::endl
    //<< dataPoints << std::endl
    //<< "Query points:" << std::endl
    //<< queryPoints << std::endl
    //<< "Neighbor indices:" << std::endl
    //<< indices << std::endl
    //<< "Neighbor distances:" << std::endl
    //<< distances << std::endl;

  for (int i=0; i<400; i++)
  {
    size_t index =  indices(i);
    auto it = map.begin();
    std::advance(it,index);
    auto t = std::find(ret.begin(),ret.end(),it->second);
    if (t==ret.end())
      ret.push_back(it->second);
  }

  return ret;
}


inline void barrel_correction(double& x, double& y, double& k1, double& k2)
{

  double rad2 = x*x + y*y;
  double rad4 = rad2*rad2;

  x*=(1.+k1*rad2 + k2*rad4);
  y*=(1.+k1*rad2 + k2*rad4);

} 

inline arma::Mat<double> create_view( double xp, double yp, double dp, double k1, double k2, arma::Mat<double>& cm, 
    lon_lat_coords& lons, lon_lat_coords& lats, country_coords& cc )
{
  arma::Mat<double> position(4,1); position.ones();
  double x,y;

  bool correct =  (fabs(k1) > 1e-6 || fabs(k2) > 1e-6) ? true : false;

  double xlimit = xp/2.;
  double ylimit = yp/2.;

  std::map<int,std::vector<lonLat>> latitudes,longitudes,countries;

  int index=0,total_objects=0;

  for (auto &object : lats)
  {
    std::vector<lonLat> tmp;
    for (auto &point: object.second)
    {    

      position(0) = point.x;
      position(1) = point.y;
      position(2) = point.z;

      arma::Mat<double> cframe = cm*position;

      double imx = -cframe(1) / cframe(2);
      double imy = cframe(0) / cframe(2);

      if (fabs(imx)<=xlimit && fabs(imy)<=ylimit)
      {
        //if (correct) barrel_correction(imx,imy,k1,k2);
        std::cerr << imx << " " << imy << " " << xlimit << " " << ylimit << std::endl;
        point.imx=imx;
        point.imy=imy;
        total_objects++;
        tmp.push_back(point);
      }
    }

    if (!tmp.empty())
    {
      latitudes.insert(std::pair<int,std::vector<lonLat> >{index,tmp});
      index++;
    }
  }

  index=0;
  for (auto &object : lons)
  {
    std::vector<lonLat> tmp;
    for (auto &point: object.second)
    {    

      position(0) = point.x;
      position(1) = point.y;
      position(2) = point.z;

      arma::Mat<double> cframe = cm*position;

      double imx = -cframe(1) / cframe(2);
      double imy = cframe(0) / cframe(2);

      if (fabs(imx)<=xlimit && fabs(imy)<=ylimit)
      {
        if (correct) barrel_correction(imx,imy,k1,k2);
        //std::cerr << imx << " " << imy << std::endl;
        point.imx=imx;
        point.imy=imy;
        total_objects++;
        tmp.push_back(point);
      }
    }

    if (!tmp.empty())
    {
      longitudes.insert(std::pair<int,std::vector<lonLat> >{index,tmp});
      index++;
    }
  }

  index = 0;
  for (auto &object : cc)
  {
    for (auto &points: object.second)
    { 
      std::vector<lonLat> tmp;
      for (auto &point: points)
      {    
        position(0) = point.x;
        position(1) = point.y;
        position(2) = point.z;

        arma::Mat<double> cframe = cm*position;

        double imx = -cframe(1) / cframe(2);
        double imy = cframe(0) / cframe(2);

        if (fabs(imx)<=xlimit && fabs(imy)<=ylimit)
        {
          if (correct) barrel_correction(imx,imy,k1,k2);
          point.imx=imx;
          point.imy=imy;
          total_objects++;
          tmp.push_back(point);
        }
      }

      if (!tmp.empty())
      {
        countries.insert(std::pair<int,std::vector<lonLat>>{index,tmp});
        index++;
      }
    }
  }

  arma::Mat<double> ret(total_objects,9);


  index = 0;
  for (auto it = longitudes.begin(); it!=longitudes.end(); it++)
  {
    int subtype = it->first; 
    for (auto &p : it->second)
    {
      ret(index,0) = 0.;
      ret(index,1) = (double) subtype;
      ret(index,2) = p.lon;
      ret(index,3) = p.lat;
      ret(index,4) = p.x;
      ret(index,5) = p.y;
      ret(index,6) = p.z;
      ret(index,7) = p.imx;
      ret(index,8) = p.imy;

      index++;
    }
  }

  for (auto it = latitudes.begin(); it!=latitudes.end(); it++)
  {
    int subtype = it->first; 
    for (auto &p : it->second)
    {
      ret(index,0) = 1.;
      ret(index,1) = (double) subtype;
      ret(index,2) = p.lon;
      ret(index,3) = p.lat;
      ret(index,4) = p.x;
      ret(index,5) = p.y;
      ret(index,6) = p.z;
      ret(index,7) = p.imx;
      ret(index,8) = p.imy;

      index++;
    }
  }

  for (auto it = countries.begin(); it!=countries.end(); it++)
  {
    int subtype = it->first; 
    for (auto &p : it->second)
    {
      ret(index,0) = 2.;
      ret(index,1) = (double) subtype;
      ret(index,2) = p.lon;
      ret(index,3) = p.lat;
      ret(index,4) = p.x;
      ret(index,5) = p.y;
      ret(index,6) = p.z;
      ret(index,7) = p.imx;
      ret(index,8) = p.imy;

      index++;
    }
  }

  return ret;

}

using namespace boost::program_options;

int main(int argc, char * argv[])
{

  double eln = 142.32,
         elt = -10.60,
         dr=1.,d=10000.,
         r=800,
         dp=4,
         xp=1000.,yp=1000.;

  double yaw=8.,pitch=0.,roll=0.;
  double k1=0.,k2=0.;
  std::string world_polys,out;
  bool show_countries=true;

  options_description desc("satCamView - create satellite camera lon/lat grid using intrinsic rotations and position information. Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("wp", value<std::string>(), "world polygons kml file")
    ("eln", value<double>(), "satellite @ earth surface, longitude degrees (default = 142.32)")
    ("elt", value<double>(), "satellite @ earth surface, latitude degrees (default = -10.60)")
    ("dr", value<double>(), "grid resolution degrees (default = 1)")
    ("a", value<double>(), "satellite yaw deg.(default = 8)")
    ("b", value<double>(), "satellite pitch deg.(default = 0)")
    ("g", value<double>(), "satellite roll deg.(default = 0)")
    ("r", value<double>(), "satellite altitude km (default = 800km)")
    ("dp", value<double>(), "distance per pixel (default = 4km)")
    ("xp", value<double>(), "image pixels in x (default=1000")
    ("yp", value<double>(), "image pixels in y (default=1000)")
    ("k1", value<double>(), "barrel correction coeff. k1 (default=0)")
    ("k2", value<double>(), "barrel correction coeff. k2 (default=0)")
    ("out", value<std::string>(), "output CSV filename for grid coordinates(else write stdout)")
    ("cov", value<std::string>(), "output grid coverage image matrix (default=cov.h5)")
    ("d", value<double>(), "Euclidean distance for kdtree (nn points search) in km (default = 10000 km )");

  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  if (vm.count("help"))
  {
    std::cerr << desc << std::endl;
    return 0;
  }

  if (vm.count("wp"))
    world_polys = vm["wp"].as<std::string>();

  else
  {
    std::cerr << "no world polygons kml filename, will just output grid" << std::endl;
    show_countries=false;
  }

  if (vm.count("elt"))
    elt = vm["elt"].as<double>();

  if (vm.count("eln"))
    eln = vm["eln"].as<double>();

  if (vm.count("d"))
    d = vm["d"].as<double>();

  if (vm.count("dr"))
    dr = vm["dr"].as<double>();

  if (vm.count("dp"))
    dp = vm["dp"].as<double>();

  if (vm.count("out"))
    out = vm["out"].as<std::string>();

  d*=1000.;
  double d2 = d*d;
  coords_country coords_country_map;
  country_coords country_coords_map;
  country_coords country_coords_map_selected;

  if (show_countries)
  {
    int result = create_luts(world_polys,coords_country_map,country_coords_map);

    if (!result) return 1;

    auto cnn = find_nns(coords_country_map, elt, eln, d2);

    if (!cnn.empty())
    {
      std::cerr << "found nn countries: " << std::endl;
      for (auto &x : cnn)
      {
        std::cerr << x << std::endl;
        auto it = country_coords_map.find(x);
        country_coords_map_selected.insert(*it);
      }
    }
  }
  
  lon_lat_coords lats, lons;
  gen_grid(elt, eln, d, dr, lons,lats);
  arma::Mat<double> cm = camera_matrix(eln, elt, r, dp, yaw, pitch, roll);
  arma::Mat<double> output_grid = create_view( xp,yp,dp,k1,k2,cm, lons,lats,country_coords_map_selected);
  
  if (!out.empty())
    output_grid.save(out.c_str(),arma::csv_ascii);


  if (vm.count("cov"))
  {
    std::string cov = vm["cov"].as<std::string>();
    
    int minx = static_cast<int>(output_grid.col(7).min());   
    int maxx = static_cast<int>(output_grid.col(7).max());   
    int miny = static_cast<int>(output_grid.col(8).min());   
    int maxy = static_cast<int>(output_grid.col(8).max());

    
    int szx = maxx-minx+2;
    int szy = maxy-miny+2;
 
    //std::cerr << szx << " " << szy << std::endl;  
    arma::Mat<double> cov_im(szx,szy,arma::fill::ones);
  
    cov_im*=255;

    for (int i=0; i<output_grid.n_rows; i++)
    {
      int ii = static_cast<int>(output_grid(i,7))-minx;
      int jj = static_cast<int>(output_grid(i,8))-miny;

      //std::cerr << ii << " " << jj << " " << szx << " " << szy << std::endl;

      cov_im(ii,jj)=0.;
      
      if (ii-1>=0)
      cov_im(ii-1,jj)=0.;
      if (ii-1>=0 && jj-1>=0)
      cov_im(ii-1,jj-1)=0.;
      if (jj-1>=0)
      cov_im(ii,jj-1)=0.;
      if (jj-1>=0 && ii+1<szx)
      cov_im(ii+1,jj-1)=0.;
      if (ii+1<szx)
      cov_im(ii+1,jj)=0.;
      if (ii+1<szx && jj+1<szy)
      cov_im(ii+1,jj+1)=0.;
      if (jj+1<szy)
      cov_im(ii,jj+1)=0.;
      if (ii-1>=0 && jj+1<szy)
      cov_im(ii-1,jj+1)=0.;

    }

      //arma::inplace_trans(cov_im);
      //arma::inplace_flipud(cov_im);
      cv::Mat out_im(szx,szy,CV_8UC1);

      for (int i=0; i<szx; i++)
        for (int j=0; j<szy; j++)
          out_im.at<unsigned char>(i,j) = static_cast<unsigned char>(cov_im(i,j));


      cv::imwrite(cov.c_str(),out_im);


    //std::cerr << minx << " " << maxx << " " << miny << " " << maxy << std::endl;
  }

  //std::cerr << output_grid << std::endl;
 arma::inplace_trans(output_grid);
 
 
  if (out.empty())
    std::cout.write((char*)output_grid.memptr(),output_grid.n_elem*sizeof(double));
  return 0;
}

