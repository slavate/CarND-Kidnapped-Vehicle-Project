/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  default_random_engine gen;
  // Set the number of particles
  num_particles = 50;
  
  // Normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  
  for (int i = 0; i < num_particles; ++i) {
    // Initialize all weights to 1.
    weights.push_back(1.0);
    
    // Initialize all particles.
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  // if I put random engine into for loop, it does not works corectly!!!
  default_random_engine gen;
  
  bool use_yaw_rate = true;
  double c1 = 1.0;
  double c2 = 1.0;
  //check division by zero
  if (fabs(yaw_rate) > 1E-6) {
    // 1st constant value
    c1 = velocity / yaw_rate;
    // 2nd constant value
    c2 = delta_t * yaw_rate;
  }
  else {
    cout << "prediction() - Warning - zero yaw rate!" << endl;
    use_yaw_rate = false;
    c1 = velocity*delta_t;
  }

  for (int i = 0; i < num_particles; ++i) {
    Particle p = particles[i];

    if (use_yaw_rate) {
      p.x = p.x + c1*(sin(p.theta + c2) - sin(p.theta));
      p.y = p.y + c1*(cos(p.theta) - cos(p.theta + c2));
      p.theta = p.theta + c2;
    }
    else {
      p.x = p.x + c1*cos(p.theta);
      p.y = p.y + c1*sin(p.theta);
    }
    // Normal (Gaussian) distribution for x, y and theta
    normal_distribution<double> dist_x(p.x, std_pos[0]);
    normal_distribution<double> dist_y(p.y, std_pos[1]);
    normal_distribution<double> dist_theta(p.theta, std_pos[2]);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

    //cout << "x " << p.x << ", y " << p.y << ", yaw " << p.theta << endl;
    //cout << "x_norm " << particles[i].x << ", y_norm " << particles[i].y << ", yaw_norm " << particles[i].theta << endl;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  
  // Transformation
  for (int i = 0; i < num_particles; ++i) {

    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
    particles[i].associations.clear();

    for (int j = 0; j < observations.size(); ++j){
      double xt = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
      double yt = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);

      particles[i].sense_x.push_back(xt);
      particles[i].sense_y.push_back(yt);
      // initialize associations
      particles[i].associations.push_back(0);
    }
  }

  // Consider only landmarks in predefined range
  std::vector<LandmarkObs> landmarks_in_range;
  for (int i = 0; i < map_landmarks.landmark_list.size(); ++i) {
    // map signals from single_landmark_s (map) structure to LandmarkObs (helpers) structure
    LandmarkObs l;
    l.id = map_landmarks.landmark_list[i].id_i;
    l.x = map_landmarks.landmark_list[i].x_f;
    l.y = map_landmarks.landmark_list[i].y_f;
    double dist = sqrt(pow(l.x, 2) + pow(l.y, 2));
    if (dist <= sensor_range) {
      landmarks_in_range.push_back(l);
    }
  }
  
  double sigma_x = std_landmark[0];
  double sigma_y = std_landmark[1];

  // Associations
  for (int i = 0; i < num_particles; ++i) {
    // Observations Multivariate-Gaussian Probability
    std::vector<double> obs_weights;
    // Go through all predicted measurements for particles
    for (int l = 0; l < particles[i].sense_x.size(); ++l) {
      // Go through all landmarks
      // temp vector with distances from particle to the landmarks
      std::vector<double> r;
      
      double x = particles[i].sense_x[l];
      double y = particles[i].sense_y[l];
      
      for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
        // Calculate euclidean distance between each landmark and predicted particle position
        r.push_back(dist(x, y,
                         map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f));
      }
      // find the best match with shortest distance
      int match_index = distance(begin(r), min_element(begin(r), end(r)));
      particles[i].associations[l] = map_landmarks.landmark_list[match_index].id_i;
      //cout << "association for particle Nr " << i << ", measurement Nr " << l << " is landmark: " << particles[i].associations[l] << endl;
      
      double mu_x = map_landmarks.landmark_list[match_index].x_f;
      double mu_y = map_landmarks.landmark_list[match_index].y_f;      
      double c1 = 1 / (2 * M_PI * sigma_x * sigma_y);
      // calculate weight
      double w = c1 * exp(-0.5 * (pow((x - mu_x), 2) / pow(sigma_x, 2) +
                                  pow((y - mu_y), 2) / pow(sigma_y, 2)));

      //avoid very small weights
      if (w < 1E-100) {
        cout << "updateWeights() - Warning - w is very small!!!" << endl;
        w = 1E-50;
      }
      
      obs_weights.push_back(w);
    }

    // Update weights
    double new_weight = 1.0;
    for (int k = 0; k < obs_weights.size(); ++k){
      new_weight *= obs_weights[k];

      if (new_weight < 1E-100) {
        cout << "updateWeights() - Warning - very small weight!!!" << endl;
      }
    }

    weights[i] = new_weight;
    particles[i].weight = new_weight;

    //cout << "particle Nr " << i << ", weight " << weights[i] << endl;
  }

  //// Normalize weights
  //double weight_sum = 0.0;
  //for (int i = 0; i < num_particles; ++i) {
  //  weight_sum += particles[i].weight;
  //}

  //for (int i = 0; i < num_particles; ++i) {
  //  particles[i].weight = particles[i].weight / weight_sum;
  //  weights[i] = particles[i].weight;
  //}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  default_random_engine gen;
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::vector<Particle> new_particles;
  for (int n = 0; n<num_particles; ++n) {
    new_particles.push_back(particles[d(gen)]);
  }
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
