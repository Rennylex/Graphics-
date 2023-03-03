#version 330 core

uniform vec2 iResolution;
uniform float iTime;
uniform int iFrame;
in vec2 fragCoord; 
out vec4 fragColor;

uniform sampler2D bufferTexture;

#define M_PI 3.1415925585

////data structures for ray tracing
struct camera{
    vec3 origin;
    vec3 horizontal;
    vec3 vertical;
    vec3 LowerLeftCorner;
};

struct ray{
    vec3 ori;
    vec3 dir;
};

struct sphere{
    vec3 ori;			////sphere center
    float r;			////sphere radius
    vec3 color;			////sphere color
};

struct light {
   vec3 position;		////point light position
   vec3 color;			////point light color
};
    
struct hit{
    float t;			////parameter in the ray function
    vec3 p;				////intersection point
    vec3 normal;		////normal on the intersection point
    vec3 color;			////color of the intersecting object
};

//////////// Random functions ///////////
float g_seed = 0.;

uint base_hash(uvec2 p) {
    p = 1103515245U*((p >> 1U)^(p.yx));
    uint h32 = 1103515245U*((p.x)^(p.y>>3U));
    return h32^(h32 >> 16);
}

void init_rand(in vec2 frag_coord, in float time) {
    g_seed = float(base_hash(floatBitsToUint(frag_coord)))/float(0xffffffffU)+time;
}

vec2 rand2(inout float seed) {
    uint n = base_hash(floatBitsToUint(vec2(seed+=.1,seed+=.1)));
    uvec2 rz = uvec2(n, n*48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU))/float(0x7fffffff);
}
/////////////////////////////////////////

const float minT = 0.001;
const float maxT = 1e8;
const int numberOfSampling = 50;

////if no hit is detected, return dummyHit
const hit dummyHit = hit(-1.0, vec3(0), vec3(0), vec3(0));

////calculate the ray for a given uv coordinate
ray getRay(camera c, vec2 uv)
{
    return ray(c.origin, c.LowerLeftCorner + uv.x * c.horizontal + uv.y * c.vertical - c.origin);
}

////TODO: implement ray-sphere intersection
hit hitSphere(const ray r, const sphere s) {
    float delta = 0.f;
	////TODO: check whether r is interescting with s by updating delta
	/*Your implementation*/

    vec3 origin = r.ori;
    vec3 direction = r.dir;
    vec3 center = s.ori;
    float radius = s.r;
    float a = dot(direction,direction);
    float b = 2 * dot(origin-center,direction);
    float c = dot(origin-center,origin-center) - (radius * radius);

    delta = b*b - 4*a*c;
    
    if(delta<0.0){
        // no solution, return dummy
        return  dummyHit;
    }

    else {
		hit h;
		h.color=s.color;

		////TODO: update other attributes of hit when an intersection is detected
		/*Your implementation*/
        float t1 = (-b + sqrt(delta)) / (2*a);
        float t2 = (-b - sqrt(delta)) / (2*a);
        
        //float t1 = (-dot(direction,origin-center)+sqrt(delta))/dot(direction,direction);
        //float t2 = (-dot(direction,origin-center)-sqrt(delta))/dot(direction,direction);

        if(t1>0.0 && t2>0.0){
            h.t = min(t1,t2);
        }
        else if(t1>0.0){
            h.t = t1;
        }
        else if(t2>0.0){
            h.t = t2;
        }
        else{
            return dummyHit;
        }

        h.p = origin + h.t * direction;
        h.normal = normalize(h.p - center);
        //make all of its element possitive


        return h;
    }
}
////TODO: return the hit sphere with the smallest t
hit findHit(ray r, sphere[4] s) 
{
	hit h = dummyHit;
    ////TODO: traverse all the spheres and find the intersecting one with the smallest t value
	/*Your implementation*/
    h.t=1e8;
    for (int i=0; i<4; i++) {
        hit temp = hitSphere(r, s[i]);
        if (temp.t > 0. && temp.t < h.t) {
            h = temp;
        }
    }
	return h;
}

////TODO: calculate the pixel color for each ray
vec3 color(ray r, sphere[4] s, light[2] l)
{
    vec3 col = vec3(0);
    hit h = findHit(r, s);
    if(h.t > 0.){
		////TODO: traverse all the lights and calculate the color contribution from each of them
		////TODO: send an additional shadow ray for each light to implement the shadow effect
		/*Your implementation*/
        // for each light
            // shoot shadow ray to detect if light not obstructed by some object (hitSphere)
            // if light not occluded, then use simple Lambertain shading to calculate its color
            // else: skip this light
            // final color of the ray is sum from all lights
        for (int i=0; i<2; i++) {
            // make shadow ray to shoot towards light source
            ray shadow_ray; //= (h.p - h.normal);
            bool isOccluded = false;
            float epsilon=0.1;
            shadow_ray.dir = normalize(l[i].position - h.p);
            shadow_ray.ori = h.p + epsilon * shadow_ray.dir;//by adding epsilon, we make sure the shadow ray is not intersecting with the object itself
            //for  each  light  you  shoot  a  shadow  ray  to  detect 
            // whether  the  light  is  occluded  by  some  other  object  in  the  scene  (here  you  may  use  the  function 
            // hitSphere again)
            for (int j=0; j<4; j++) {
                hit temp = hitSphere(shadow_ray, s[j]);
                if (temp.t > 0.0 ) {
                    isOccluded = true;
                    j=4;
                }
            }


            // check if light is not occluded
            if (isOccluded == false) {
                    // ////Lambertian lighting calculation
                    // /* variables */
                    vec3 final = vec3(1);
                    vec3 light_color=l[i].color;

	                // float ka = 0.1;
	                // float kd = 0.;
	                // // declarations
	                // vec3 ambient;
	                // vec3 diffuse;
	                // ////Lambertian = Ambient + Diffuse
	                // // ambient
	                // ambient = ka * Ia;
	                // // diffuse
	                // vec3 l_dir = normalize(l[i].position - h.p);
	                // diffuse = kd * Id * max(dot(h.normal,l_dir), 0.0);
	                // final = (ambient + diffuse)*h.color;
	                // //col += final;
                    // //col=h.normal;


                    	vec3 LightPosition = l[i].position;
                        const vec3 LightIntensity = vec3(8);
                        const vec3 ka = 0.8*vec3(1., 1., 1.);
                        const vec3 kd = 0.8*vec3(1., 1., 1.);
                        const vec3 ks = vec3(2.);
                        const float n = 400.0;

                        vec3 normal = h.normal;
                        vec3 lightDir = l[i].position - h.p;
                        float _lightDist = length(lightDir);
                        vec3 _lightDir = normalize(lightDir);
                        vec3 _localLight = LightIntensity / (_lightDist * _lightDist);
                        vec3 ambColor = ka;
                        vec3 lightingColor = kd * _localLight * max(0., dot(_lightDir, normal));
                        final=h.color*lightingColor;
                    col+=final;

                }
        }
            //       In color, you will calculate the ray color for a given hit by looking at the existing light sources in
            // the  scene.  The  logic  of  this  function  is  simple:  for  each  light  you  shoot  a  shadow  ray  to  detect
            // whether  the  light  is  occluded  by  some  other  object  in  the  scene  (here  you  may  use  the  function
            // hitSphere again). If the light is not occluded, then you use a simple Lambertian shading model to
            // calculate its color; otherwise you skip this light. The final color of the ray is the sum from all the
            // lights (the number is hard-coded as 2 here).
    }
    //return col;

    return col;

}



vec3 color_(ray r, sphere[4] s, light[2] l)
{
    vec3 col = vec3(0);
    hit h = findHit(r, s);
    vec3 final = vec3(1);
    if(h.t > 0.){
		////TODO: traverse all the lights and calculate the color contribution from each of them
		////TODO: send an additional shadow ray for each light to implement the shadow effect
		/*Your implementation*/
        // for each light
            // shoot shadow ray to detect if light not obstructed by some object (hitSphere)
            // if light not occluded, then use simple Lambertain shading to calculate its color
            // else: skip this light
            // final color of the ray is sum from all lights
        for (int i=0; i<2; i++) {
            // make shadow ray to shoot towards light source
            ray shadow_ray; //= (h.p - h.normal);
            shadow_ray.ori = l[i].position;
            shadow_ray.dir = h.p - l[i].position+0.2;
            for (int j=0; j<4; j++) {
                hit h2 = hitSphere(shadow_ray, s[j]);
                // check if light is not occluded
                if (h2 == dummyHit) {
                    ////Lambertian lighting calculation
                    /* variables */
                	vec3 Ia = vec3(1, 1, 1);
	                vec3 Id = vec3(1, 1, 1);
	                float ka = 0.2;
	                float kd = 0.9;
	                // declarations
	                vec3 ambient;
	                vec3 diffuse;
	                ////Lambertian = Ambient + Diffuse
	                // ambient
	                ambient = ka * Ia;
	                // diffuse
	                vec3 l_dir = normalize(l[i].position - h.p);
	                diffuse = kd * Id * max(dot(h.normal, l_dir), 0.0);
	                final = ambient + diffuse;
                    
	                //col += final * h.color;
                }
                else
                    final = vec3(1);
                    
            }
            col += final * h.color;
        }
 //       In color, you will calculate the ray color for a given hit by looking at the existing light sources in
 // the  scene.  The  logic  of  this  function  is  simple:  for  each  light  you  shoot  a  shadow  ray  to  detect
// whether  the  light  is  occluded  by  some  other  object  in  the  scene  (here  you  may  use  the  function
// hitSphere again). If the light is not occluded, then you use a simple Lambertian shading model to
// calculate its color; otherwise you skip this light. The final color of the ray is the sum from all the
// lights (the number is hard-coded as 2 here).
    }
    //return col;

    return col;

}





void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    camera c = camera(vec3(0,15,50), vec3(5, 0, 0), vec3(0, 3, -3), vec3(-2.5, -1.5, -1));
    sphere s[4];
    s[0] = sphere(vec3(0, 0.6, -1), 0.6, vec3(0.8,0.2,0.2));
	s[1] = sphere(vec3(1.2, 0.4, -1), 0.4, vec3(0.2,0.9,0.2));
	s[2] = sphere(vec3(-1.2, 0.5, -1), 0.5, vec3(0.2,0.2,0.9));
    s[3] = sphere(vec3(0, -200, -1),200.0, vec3(0.5,0.5,0.5));

	light l[2];
	l[0] = light(vec3(-1, 3, 0.5), vec3(1));
	l[1] = light(vec3(0.5, 2, 1), vec3(1));
    vec3 resultCol = vec3(0);

    // Here I use i to get differnet seeds for each run
    init_rand(fragCoord, iTime);
    vec2 random = rand2(g_seed);
    ray r = getRay(c, uv + random/iResolution.xy);
    resultCol += color(r, s, l);
    
	// Output to screen
    fragColor = vec4((resultCol + float(iFrame-1) * texture(bufferTexture, uv).xyz)/float(iFrame), 1.);
}

void main() {
	mainImage(fragColor, fragCoord);
}
