#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement MirrorBSDF
  *pdf = 1.0;
  reflect(wo, wi);
  return reflectance / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Assignment 7: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  
  double roughness = alpha,
         numer = exp(-pow(tan(getTheta(h)) / roughness, 2)),
         denom = PI * pow(roughness, 2) * pow(cos(getTheta(h)), 4);

  return numer / denom;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
  double cos_theta_i = cos(getTheta(wi));
  Vector3D Rs = ((eta*eta + k*k) - 2*eta*cos_theta_i + pow(cos_theta_i, 2)) /
                ((eta*eta + k*k) + 2*eta*cos_theta_i + pow(cos_theta_i, 2)),
           Rp = ((eta*eta + k*k) * pow(cos_theta_i, 2) - 2*eta*cos_theta_i + 1) /
                ((eta*eta + k*k) * pow(cos_theta_i, 2) + 2*eta*cos_theta_i + 1); 
  return (Rs + Rp) / 2;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Implement microfacet model here.
  Vector3D h = (wi + wo).unit();
  
  if (wo.z > 0 && wi.z > 0)
    return F(wi) * G(wo, wi) * D(h) / (4 * wo.z * wi.z);

  return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  Vector2D rand = sampler.get_sample();
  double r1 = rand.x,
         r2 = rand.y,
         theta_h = atan(sqrt(-pow(alpha, 2) * log(1-r1))),
         phi_h = 2 * PI * r2;
  
  Vector3D h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
  *wi = wo + 2 * (dot(wo, h)*h - wo);

  double p_theta = 2 * sin(theta_h) * exp(-pow(tan(theta_h), 2) / pow(alpha, 2)) / 
                   (pow(alpha, 2) * pow(cos(theta_h), 3)),
         p_phi = 1.0 / (2 * PI),
         p_w_h = p_theta * p_phi / sin(theta_h),
         p_w_wi = p_w_h / (4 * dot((*wi), h));

  if (wi->z > 0 && wo.z > 0) {
    *pdf = p_w_wi;
    return MicrofacetBSDF::f(wo, *wi);
  }
  *wi = cosineHemisphereSampler.get_sample(pdf);
  return Vector3D();
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement RefractionBSDF

  if (refract(wo, wi, ior)) {
    *pdf = 1.0;
    double eta = wo.z > 0 ? 1.0 / ior : ior;
    return transmittance / abs_cos_theta(*wi) / (eta*eta);
  }
  return Vector3D();
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

  double eta = wo.z > 0 ? 1.0 / ior : ior,
         check_total_ref = 1 - eta*eta*(1-wo.z*wo.z);
  
  if (check_total_ref < 0) {
    *wi = wo;
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  }

  refract(wo, wi, ior); // give value to wi
  
  double R0 = pow((1-eta) / (1+eta), 2),
         R = R0 + (1-R0) * pow(1 - abs_cos_theta(*wi), 5);

  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return R * reflectance / abs_cos_theta(*wi);
  }
  else {
    refract(wo, wi, ior);
    *pdf = 1 - R;
    return (1-R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
  }
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {
  // TODO Assignment 7: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.

  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Assignment 7: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  // entering non-air object
  double eta = wo.z > 0 ? 1.0 / ior : ior,
           k = 1 - eta*eta*(1-wo.z*wo.z);
  if (k < 0)
    return false;

  *wi = Vector3D(-eta * wo.x, -eta * wo.y, (wo.z < 0 ? 1.0 : -1.0) * sqrt(k));

  return true;
}

} // namespace CGL
