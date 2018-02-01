#version 300 es

precision highp float;

uniform vec2 u_Dimensions;
uniform float u_Fov;
uniform float u_Time;
uniform vec3 u_EyePos;
uniform mat4 u_ViewInv; 
uniform sampler2D u_BackgroundTexture;

//----------------------------for background------------------
const vec3 flame[5] = vec3[](vec3(255, 102, 0) / 255.0,
							vec3(255, 154, 0) / 255.0,
							vec3(0.0, 0.0, 0.0),
							vec3(0.0, 0.0, 0.0),
							vec3(0.0, 0.0, 0.0));

const float PI = 3.14;

vec2 sphereToUV(vec3);
vec3 uvToDusk(vec2);

float fbm(const in vec2 uv);
float noise(in vec2 uv);
vec2 smoothF(vec2 uv);
//----------------------------------------------------------------------

out vec4 out_Col;

mat4 translationMatrix(vec3 v) {
	return mat4(1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                v.x, v.y, v.z, 1.0);		
}

mat4 rotateX(float rad) {
	return mat4(1.0, 0.0, 0.0, 0.0,
                0.0, cos(rad), sin(rad), 0.0,
                0.0, -sin(rad), cos(rad), 0.0,
				0.0, 0.0, 0.0, 1.0);
}

mat4 rotateY(float rad) {
	return mat4(cos(rad), 0.0, -sin(rad), 0.0,
                0.0, 1.0, 0.0, 0.0,
                sin(rad), 0.0, cos(rad), 0.0,
				0.0, 0.0, 0.0, 1.0);
}

mat4 rotateZ(float rad) {
	return mat4(cos(rad), sin(rad), 0.0, 0.0,
                -sin(rad), cos(rad), 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0);
}


float length2(vec2 p) {
	return pow((pow(p.x, 2.0) + pow(p.y, 2.0)), 1.0 / 2.0);
}

float length8(vec2 p) {
	return pow((pow(p.x, 8.0) + pow(p.y, 8.0)), 1.0 / 8.0);
}

// -------------------------SDF operations--------------------
float smin(float a, float b, float k) {
	float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
	return mix(b, a, h) - k * h * (1.0 - h);
}

float unionSDF(float a, float b) {
	return min(a, b);
}


// -------------------------primitive SDF--------------------
float sdTorus82( vec3 p, vec2 t )
{
  vec2 q = vec2(length2(p.xz)-t.x,p.y);
  return length8(q)-t.y;
}

float sdPlane( vec3 p, vec4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float udRoundBox( vec3 p, vec3 b, float r )
{
  return length(max(abs(p)-b,0.0))-r;
}

float sdCone( vec3 p, vec2 c )
{
    // c must be normalized
    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}

float charO(vec3 p, vec2 t) {
  mat4 transM = rotateX(radians(u_Time)) * rotateY(radians(u_Time)) * rotateY(radians(u_Time));
  vec3 transedP =  vec3(transM * vec4(p, 1.0));
  vec2 q = vec2(length2(transedP.xz)-t.x,transedP.y);
  
  return length8(q)-t.y;
}

float charE(vec3 p) {
	mat4 transM = rotateX(radians(u_Time)) * rotateY(radians(u_Time)) * rotateY(radians(u_Time));
  	vec3 transedP =  vec3(transM * vec4(p, 1.0));
   
    transM = translationMatrix(vec3(-3.0, 0.0, 0.0));
	float left = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(1.0, 5.0, 1.0), 0.5);

	transM = translationMatrix(vec3(0.0, 4.5, 0.0));
	float up = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(4.0, 1.0, 1.0), 0.5);
	
	transM = translationMatrix(vec3(0.0, 0.0, 0.0));
	float middle = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(4.0, 1.0, 1.0), 0.5); 

	transM = translationMatrix(vec3(0.0, -4.5, 0.0));
	float bottom = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(4.0, 1.0, 1.0), 0.5);

	return unionSDF(unionSDF(unionSDF(up, middle), bottom), left);
}

float charH(vec3 p) {
	mat4 transM = rotateX(radians(u_Time)) * rotateZ(radians(u_Time)) * rotateY(radians(u_Time));
	vec3 transedP =  vec3(transM * vec4(p, 1.0));
	transM = translationMatrix(vec3(-3.0, 0.0, 0.0));
	float left = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(1.0, 5.0, 1.0), 0.5);

	transM = translationMatrix(vec3(3.0, 0.0, 0.0));
	float right = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(1.0, 5.0, 1.0), 0.5);

	transM = translationMatrix(vec3(0.0, 0.0, 0.0));
	float middle = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(4.0, 1.0, 1.0), 0.5);

	return unionSDF(unionSDF(left, right), middle);
}

float charL(vec3 p) {
	mat4 transM = rotateY(radians(u_Time)) * rotateX(radians(u_Time)) * rotateZ(radians(u_Time));
	vec3 transedP =  vec3(transM * vec4(p, 1.0));

	transM = translationMatrix(vec3(-3.0, 0.0, 0.0));
	float vert = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(1.0, 5.0, 1.0), 0.5);

	transM = translationMatrix(vec3(0.0, -4, 0.0));
	float hori = udRoundBox(vec3(inverse(transM) * vec4(transedP, 1.0)), vec3(4.0, 1.0, 1.0), 0.5);
	return unionSDF(vert, hori);
}

float mirror(vec3 p) {
	float mirror = udRoundBox(p, vec3(15.0, 0.2, 15.0), 0.5);
	return mirror;
}

float backMirror(vec3 p) {
	
	float mirror = udRoundBox(p, vec3(20.0, 20.0, 1.0), 0.5);

	mat4 transM = rotateX(radians(90.0));
	vec3 transedP =  vec3(transM * vec4(p, 1.0));
	float torus = sdTorus82(transedP, vec2(7.0, 7.0));

	vec3 c = vec3(2.0, 0.0, 2.0);
	transedP = mod(transedP, c) - 0.5 *c;
	float ball = sdSphere(transedP, 1.0);

	return max(-max(-torus, mirror), ball);
}


// -------------------------Scene SDF--------------------
float sceneSDF(vec3 samplePoint) {
	float[7] SDFArray;
	
	// H char
	mat4 transM = translationMatrix(vec3(20.0 * sin(0.05 * u_Time) - 20.0, -5.0, 0.0));
	vec3 transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	float h = charH(transedP);
	SDFArray[0] = h;

	// E char
	transM = translationMatrix(vec3(10.0 * sin(0.05 * u_Time) - 10.0, -5.0, 0.0));
	transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	float e = charE(transedP);
	SDFArray[1] = e;

	// L char
	transM = translationMatrix(vec3(0.0, -5.0, 0.0));
	transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	float l_1 = charL(transedP);
	SDFArray[2] = l_1;

	// L char
	transM = translationMatrix(vec3(10.0 * cos(0.05 * u_Time + 3.14 / 2.0) + 10.0, -5.0, 0.0));
	transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	float l_2 = charL(transedP);
	SDFArray[3] = l_2;

	// O char
	transM = translationMatrix(vec3(20.0 * cos(0.05 * u_Time + 3.14 / 2.0) + 20.0, -5.0, 0.0));
	transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	vec2 sdTorus82_Para = vec2(4.0, 1.5);
	float o = charO(transedP, sdTorus82_Para);
	SDFArray[4] = o;

	// mirror
	transM = translationMatrix(vec3(0.0, 5.0, -30.0));
	transedP = vec3(inverse(transM) * vec4(samplePoint, 1.0));  // rotation and translation
	float backMirror = backMirror(transedP);
	SDFArray[5] = backMirror;


	//plane--------------
	transM = translationMatrix(vec3(0.0, -25.0, 0.0));	
	float plane =  sdPlane(vec3(inverse(transM) * vec4(samplePoint, 1.0)), vec4(0.0, 1.0, 0.0, 1.0));
	SDFArray[6] = plane;

	float result = SDFArray[0];
	for (int i = 1; i < 7; ++i) {
		result = smin(result, SDFArray[i], 6.0);
	}
	return result;
}

vec3 rayDirection() {
    vec2 xy = vec2(gl_FragCoord.x, gl_FragCoord.y) - u_Dimensions / 2.0;
    float z = (u_Dimensions.y / 2.0) / tan(radians(u_Fov));
    return normalize(vec3(xy, -z));
}

float march(vec3 rDir, float min_t, float max_t, float EPSILON) {
	const int MAX_MARCHING_STEPS = 255;

	float t = min_t;
	for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
		vec3 p = u_EyePos + t * rDir;
    	float dist = sceneSDF(p);
    	if (dist < EPSILON) {
        	// We're inside the scene surface!
        	return t;
    	}
    	// Move along the view ray
    	t += dist;

    	if (t >= max_t) {
        // Gone too far; give up
        	return max_t;
    	}
	}
	return max_t;
}

vec3 estimateNormal(vec3 p) {
	const float EPSILON = 0.0001;
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

void main() {
	const float MIN_T = 0.0;
	const float MAX_T = 100.0;
	const float EPSILON = 0.0001;

	vec3 rDir = rayDirection();
	//vec3 rDir_transed = vec3(u_ViewInv * vec4(rDir, 1.0));
	vec3 rDir_transed = rDir;
    float t = march(rDir_transed, MIN_T, MAX_T, EPSILON);
    vec3 pos = u_EyePos + t * rDir_transed;
	vec3 normal = estimateNormal(pos);

	// lambert
	vec3 lightPos = vec3(5.0, 5.0, 3.0);
	vec3 lightVector = normalize(lightPos - pos);
	vec4 diffuseColor = vec4(1.0, 0.0, 0.0, 1.0);
    float diffuseTerm = dot(normal, lightVector);
	float ambientTerm = 0.2;
	float lightIntensity = diffuseTerm + ambientTerm; 

	if (t >= MAX_T) {
        // Didn't hit anything
        out_Col = vec4(0.0, 0.0, 0.0, 1.0);
		//-----------------------------------draw background-----------------------
		 vec2 uv = sphereToUV(rDir);

    	vec2 uvT1 = uv + vec2(u_Time * 0.01);
    	vec2 uvT2 = uv + vec2(u_Time * 0.0005, -u_Time * 0.002);

		vec2 slope = vec2(fbm(uvT2 + vec2(1.0/u_Dimensions.x, 0)) - fbm(uvT2 - vec2(1.0/u_Dimensions.x, 0)),
						fbm(uvT2 + vec2(0, 1.0/u_Dimensions.y)) - fbm(uvT2 - vec2(0, 1.0/u_Dimensions.y)));


		vec3 distortedDuskHue = uvToDusk(uv + slope);
		out_Col = vec4(distortedDuskHue, 1.0);
	//---------------------------------------------------------------------------------------------------

    } else {
		// cast another ray
		vec3 newRayDir = rDir - 2.0 * dot(rDir, normal) * normal;
		float newT = march(newRayDir, MIN_T, MAX_T, EPSILON);
		vec3 newPos = pos + newT * newRayDir;
		vec3 newNormal = estimateNormal(newPos);

		float newDiffuseTerm = dot(newNormal, lightVector);
		float newAmbientTerm = 0.2;
		float newLightIntensity = newDiffuseTerm + newAmbientTerm; 

		vec4 refColor = vec4(diffuseColor.rgb * newLightIntensity, diffuseColor.a);
		out_Col = refColor * 0.8;
	}
}

//-------------------------------background functions-----------------

vec3 uvToDusk(vec2 uv)
{
    // Below horizon
    if(uv.y < 0.5)
    {
        return flame[0];
    }
    else if(uv.y < 0.55) // 0.5 to 0.55
    {
        return mix(flame[0], flame[1], (uv.y - 0.5) / 0.05);
    }
    else if(uv.y < 0.6)// 0.55 to 0.6
    {
        return mix(flame[1], flame[2], (uv.y - 0.55) / 0.05);
    }
    else if(uv.y < 0.65) // 0.6 to 0.65
    {
        return mix(flame[2], flame[3], (uv.y - 0.6) / 0.05);
    }
    else if(uv.y < 0.75) // 0.65 to 0.75
    {
        return mix(flame[3], flame[4], (uv.y - 0.65) / 0.1);
    }
    return flame[4]; // 0.75 to 1
}

// Convert a point on a sphere to a UV coordinate
vec2 sphereToUV(vec3 p)
{
    float phi = atan(p.z, p.x); // Returns atan(z/x)
    if(phi < 0.0)
    {
        phi += 2.0 * PI; // [0, TWO_PI] range now
    }
    // ^^ Could also just add PI to phi, but this shifts where the UV loop from X = 1 to Z = -1.
    float theta = acos(p.y); // [0, PI]
    return vec2(1.0 - phi / 2.0 * PI, 1.0 - theta / PI);
}


vec2 smoothF(vec2 uv)
{
    return uv*uv*(3.-2.*uv);
}

float noise(in vec2 uv)
{
    const float k = 257.;
    vec4 l  = vec4(floor(uv),fract(uv));
    float u = l.x + l.y * k;
    vec4 v  = vec4(u, u+1.,u+k, u+k+1.);
    v       = fract(fract(1.23456789*v)*v/.987654321);
    l.zw    = smoothF(l.zw);
    l.x     = mix(v.x, v.y, l.z);
    l.y     = mix(v.z, v.w, l.z);
    return    mix(l.x, l.y, l.w);
}

float fbm(const in vec2 uv)
{
    float a = 0.5;
    float f = 5.0;
    float n = 0.;
    int it = 8;
    for(int i = 0; i < 32; i++)
    {
        if(i<it)
        {
            n += noise(uv*f)*a;
            a *= .5;
            f *= 2.;
        }
    }
    return n;
}
