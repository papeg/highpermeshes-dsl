#ifndef DSL_MESHES_DIMENSIONS_HPP
#define DSL_MESHES_DIMENSIONS_HPP

#include <HighPerMeshes/common/DataTypes.hpp>
#include <HighPerMeshes/auxiliary/Math.hpp> 

#include <map>
#include <iterator>
#include <algorithm>


namespace HPM::mesh
{
    using Vec3D = HPM::dataType::Vec3D;
   
    class DomainDimensions
    {
        public:
        Vec3D min;
        Vec3D max;
        Vec3D tfsfCenter;
        Vec3D tfsfDimension;
        Vec3D tfsfBegin;
        Vec3D tfsfEnd;
        double pmlThickness;
        Vec3D pmlBegin;
        Vec3D pmlEnd;
        
        Vec3D getCenter() const
        {
            return Vec3D{
                (max.x + min.x) / 2,
                (max.y + min.y) / 2,
                (max.z + min.z) / 2
            };
        }

        bool CoordinateIsInCenterWithRadius(Vec3D coord, double radius) const
        {
            return CoordinateIsInRadiusOfOther(coord, radius, getCenter());
        }

        bool CoordinateIsInRadiusOfOther(Vec3D coord, double radius, Vec3D other) const
        {
            return (coord.x - other.x) * (coord.x - other.x)
                + (coord.y - other.y) * (coord.y - other.y)
                + (coord.z - other.z) * (coord.z - other.z)
                <= radius * radius;    
        }

        friend std::istream& operator>>(std::istream& is, DomainDimensions &dim)
        {
            std::string line;
            while (std::getline(is, line) && (line.find("Domain dimensions:") == std::string::npos)) {}
            is >> dim.min[0] >> dim.min[1] >> dim.min[2];
            is >> dim.max[0] >> dim.max[1] >> dim.max[2];
            while (std::getline(is, line) && (line.find("TFSF parameters") == std::string::npos)) {}
            is >> dim.tfsfCenter[0] >> dim.tfsfCenter[1] >> dim.tfsfCenter[2];
            is >> dim.tfsfDimension[0] >> dim.tfsfDimension[1] >> dim.tfsfDimension[2];
            while (std::getline(is, line) && (line.find("PML thickness") == std::string::npos)) {}
            is >> dim.pmlThickness;
            for (unsigned i = 0; i < 3; i++) {
                dim.pmlBegin[i] = dim.min[i] + dim.pmlThickness;
            }
            for (unsigned i = 0; i < 3; i++) {
                dim.pmlEnd[i] = dim.max[i] - dim.pmlThickness;
            }
            for (unsigned i = 0; i < 3; i++) {
                dim.tfsfBegin[i] = dim.tfsfCenter[i] - (dim.tfsfDimension[i] / 2.0);
            }
            for (unsigned i = 0; i < 3; i++) {
                dim.tfsfEnd[i] = dim.tfsfCenter[i] + (dim.tfsfDimension[i] / 2.0);
            }
            return is;
        }

        friend std::ostream& operator<<(std::ostream& os, const DomainDimensions &dim)
        {
            os << "       MY HEADER" << std::endl;
            os << dim.min[0] << " " << dim.min[1] << " " << dim.min[2] << std::endl;
            os << dim.max[0] << " " << dim.max[1] << " " << " " << dim.max[2] << std::endl;
            os << "TFSF parameters (center and dimensions):" << std::endl;
            os << dim.tfsfCenter[0] << " " << dim.tfsfCenter[1] << " " << dim.tfsfCenter[2] << std::endl;
            os << dim.tfsfDimension[0] << " " << dim.tfsfDimension[1] << " " << dim.tfsfDimension[2] << std::endl;
            os << " PML thickness:" << std::endl;
            os << dim.pmlThickness << std::endl;
            return os;
        }
    };


}
#endif