struct vec3
{
        double x,y,z;
	vec3(double X,double Y,double Z)
	{
		x = X;
		y = Y;
		z = Z;
	}
        vec3 operator+(const vec3 other)
        {   
                return vec3(x+other.x,y+other.y,z+other.z);
        }   
        vec3 operator-(const vec3 other)
        {
                return vec3(x-other.x,y-other.y,z-other.z) ;
        }
	vec3 operator*(const int other)
	{
		return vec3(x*other,y*other,z*other);
	}

};
int CheckLineBox( vec3 B1, vec3 B2, vec3 L1, vec3 L2);
