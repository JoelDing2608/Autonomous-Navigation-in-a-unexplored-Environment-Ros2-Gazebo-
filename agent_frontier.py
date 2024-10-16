"""
Partiel 2023 Sp423
Quentin D

incertid
"""
import numpy as np
import math as m

# Parametres pour la Terre et Jupiter (à adapter selon la mission)
mu_earth = 398600  # Parametre gravitationnel de la Terre (km^3/s^2)
r_earth = 6378  # Rayon de la Terre (km)
gE = 9.81
ISP = 180

# ! orbite 1 //////////////////////////////////////////////////////////////////
print("Orbite 1 //////////////////////////////////////////////////////////")    
#Module Alt Values:
Alt_p = 480 #km
Alt_a = 650 #km

#Rayon p et a
r_p1 = Alt_p +r_earth
r_a1 = Alt_a +r_earth
print("Rayon de l'orbite r_p1 = ",r_p1,"  r_a1 = ",r_a1)

#semi major axis 
a_1 = (r_p1 + r_a1)/2
print("Demi grand axe a_1 = ",a_1)

#L’excentricité e
e_1 = (r_a1 - r_p1)/(r_a1 + r_p1)
print("excentricité e : e_1 = ",e_1)

#vitesse a p et a
v_p1 = m.sqrt(2*mu_earth*(-1/(2*a_1)+1/r_p1))
v_a1 = m.sqrt(2*mu_earth*(-1/(2*a_1)+1/r_a1))
print("Vitesse en orbite v_p1 = ",v_p1,"  v_a1 = ",v_a1)

#Time delay between launch and injection in orbit
tL = 0  #sec
delta_t0 = 980 #sec
t0 = tL + delta_t0  #sec
print("time delay to = ",t0)


# ! orbite 2 //////////////////////////////////////////////////////////////////
print("Orbite 2 //////////////////////////////////////////////////////////")
#Module Alt Values:
Alt_p2 = 650 #km
Alt_a2 = 1790 #km

#Rayon p et a
r_p2 = Alt_p2 +r_earth
r_a2 = Alt_a2 +r_earth
print("Rayon de l'orbite r_p1 = ",r_p2,"  r_a1 = ",r_a2)

#semi major axis 
a_2 = (r_p2 + r_a2)/2
print("Demi grand axe a_2 = ",a_2)

#L’excentricité e
e_2 = (r_a2 - r_p2)/(r_a2 + r_p2)
print("excentricité e : e_2 = ",e_2)

#vitesse a p et a
v_p2 = m.sqrt(2*mu_earth*(-1/(2*a_2)+1/r_p2))
v_a2 = m.sqrt(2*mu_earth*(-1/(2*a_2)+1/r_a2))
print("Vitesse en orbite v_p2 = ",v_p2,"  v_a2 = ",v_a2)

#Time delay between launch and first maneuver
delta_t1 = 23.5*3600 #sec
t1 = t0 + delta_t1  #sec
print("time delay t1 = ",t1)

################################################################################# a voir si c'est juste
# Increment of velocity avec correction a l'apoastre 
dV1_a = ((r_p2-r_p1)*mu_earth)/(4*(a_2**2)*v_a2)
#dV1_a = v_a2-v_p1
print("Increment of velocity a l'apoastre : dV1_a = ",dV1_a)
# Increment of velocity avec correction au periastre 
dV1_p = ((r_a2-r_a1)*mu_earth)/(4*(a_2**2)*v_p2)        # ? a revoir 
#dV1_p = v_a2-v_p1
print("Increment of velocity au periastre : dV1_p = ",dV1_p)

# ? avant inversion /////////////////////////////////
Alt_p2M = 650 #km
Alt_a2M = 650 #km
r_p2M = Alt_p2M +r_earth
r_a2M = Alt_a2M +r_earth
a_2M = (r_p2M + r_a2M)/2
dV1_aM = ((r_p2M-r_p1)*mu_earth)/(4*(a_2M**2)*v_a1)

# ? Apres inversion /////////////////////////////////
Alt_p2F = 650 #km
Alt_a2F = 1790 #km
r_p2F = Alt_p2F +r_earth
r_a2F = Alt_a2F +r_earth
a_2F = (r_p2F + r_a2F)/2
dV1_pF = ((r_a2F-r_a1)*mu_earth)/(4*(a_2F**2)*v_p2)
VP1_tot = dV1_aM + dV1_pF
print("Increment of velocity au periastre : VP1_tot = ",VP1_tot)
########################################################################


#ratio of masses mi/mf
ratio = 1/m.exp(dV1_a/(gE*ISP))
print("ratio of masses mi/mf : ratio = ",ratio)

# ! orbite 3 /////////////////////////////////////////////////////////////////////////
print("Orbite 3 //////////////////////////////////////////////////////////////")
#Module Alt Values:
Alt_p3 = 650 #km
Alt_a3 = 394500 #km

#Rayon p et a
r_p3 = Alt_p3 +r_earth
r_a3 = Alt_a3 +r_earth
print("Rayon de l'orbite r_p3 = ",r_p3,"  r_a3 = ",r_a3)

#semi major axis 
a_3 = (r_p3 + r_a3)/2
print("Demi grand axe a_3 = ",a_3)

#L’excentricité e
e_3 = (r_a3 - r_p3)/(r_a3 + r_p3)
print("excentricité e : e_3 = ",e_3)

#vitesse a p et a
v_p3 = m.sqrt(2*mu_earth*(-1/(2*a_3)+1/r_p3))
v_a3 = m.sqrt(2*mu_earth*(-1/(2*a_3)+1/r_a3))
print("Vitesse en orbite v_p3 = ",v_p3,"  v_a3 = ",v_a3)

#Time delay between launch and maneuve
Period = 2*m.pi*m.sqrt(a_3**3/mu_earth) # ?  
delta_t2 = Period*2 #sec
t2 = t1 + delta_t2  #sec
print("time delay t2 = ",t2)

# Increment of velocity avec correction a l'apoastre 
dV2_a = ((r_p3-r_p2)*mu_earth)/(4*(a_3**2)*v_a3)
#dV1_a = v_a2-v_p1
print("Increment of velocity a l'apoastre : dV2_a = ",dV2_a)
# Increment of velocity avec correction au periastre 
dV2_p = ((r_a3-r_a2)*mu_earth)/(4*(a_3**2)*v_p3)
#dV1_p = v_a2-v_p1
print("Increment of velocity au periastre : dV2_p = ",dV2_p)

#ratio of masses mi/mf
ratio2 = 1/m.exp(dV2_p/(gE*ISP))
print("ratio of masses mi/mf : ratio2 = ",ratio2)

# ! orbite 4 /////////////////////////////////////////////////////////////////////////
print("Orbite 4 //////////////////////////////////////////////////////////////")

mu_moon = 4903  
V_moon = 1.02 #km/s
r_Tearth = 318100
r_Tmoon = 10427
Influence_radius_moon = 66100
RM=1738 
#Velocity at the entry of the influence sphere of the Moon, from Earth ####### Verif les formules
V_entry = m.sqrt(2*mu_earth*(-1/(2*a_3)+1/(r_Tearth)))
print("Velocity at the entry of the influence sphere of the Moon, from Earth v_entry",V_entry)

gamma1=m.acos((v_p3*r_p3)/(V_entry*(r_Tearth)))
print("gamma 1=",gamma1,"rad",gamma1*180/m.pi,"°")
#velocity at the entry of the influence sphere of the Moon, from Moon ####### Verif les formules
V_inf = m.sqrt(V_entry**2+V_moon**2-2*V_entry*V_moon*m.cos(gamma1))

print("Velocity at the entry of the influence sphere of the Moon, from Earth v_entry :",V_entry)
print("Velocity at the entry of the influence sphere of the Moon, from Moon V_inf",V_inf)


#semie major axis ################## Verif les formules pas la meme que pour les autres a 
a_4 = ((Influence_radius_moon*mu_moon)/((Influence_radius_moon*V_inf**2)-(mu_moon*2)))
print("semie major axis a_4 = ",a_4)  

#Eccentricity at periselene e_4  r=a(e*cosh(E)-1) (hyperbole)

e_4 = (r_Tmoon+a_4)/a_4
print("Eccentricity at periselene e_4 = ",e_4)


phi1=m.acos(1/e_4)
print("phi 1=", phi1,"rad", phi1*180/m.pi,"°")

# velocity at the periselene V_p4
vp4=m.sqrt(2*mu_moon*(1/(2*a_4)+(1/(r_Tmoon+RM))))
print("Vp4=",vp4,"km/s")

E1 = m.acos((1/e_3)*(1-r_Tearth/a_3)) 
print("E1=", E1,"rad")
E2 = m.acosh((1/e_4)*((Influence_radius_moon/a_4)+1))
print("E2=", E2,"rad")
# Total time delay t_3 between launch and Periselene position 
# Formula  : t-tp = sqrt(a^3/mu) * (E - e sin(E))
t3_1 =  m.sqrt(a_3**3/mu_earth)*(E1 - e_3*m.sin(E1))
t3_2 =  m.sqrt(a_4**3/mu_moon)*(e_4*m.sinh(E2) - E2)

t3 = t3_1 + t3_2 + t2
print("Total time delay t_3 between launch and Periselene position t3 = ",t3,"sec")
print("t_3 in days : ",t3/60/60/24,"days")


"""
mu_jupiter = 1.267e17  # Paramètre gravitationnel de Jupiter (m^3/s^2)
r_jupiter = 69911e3  # Rayon de Jupiter (m)
# Paramètres
delta_v_initial = 7.92e3  # Delta-v initial (m/s)
a1 = 3173.26e3  # Demi-grand axe de l'hyperbole (m)

# Calcul de la vitesse de transfert (manœuvre de Hohmann)
a_transfer = (r_earth + r_jupiter) / 2  # Demi-grand axe de l'orbite de transfert
v_initial = np.sqrt(mu_earth / r_earth)  # Vitesse à l'orbite terrestre
v_final = np.sqrt(mu_jupiter / r_jupiter)  # Vitesse à l'orbite de Jupiter

# Premier delta-v
delta_v2 = np.sqrt(mu_jupiter / a_transfer) - v_final  # Deuxième delta-v

print(f"Delta-v pour la manœuvre de départ : {delta_v1:.2f} m/s")
print(f"Delta-v pour la manœuvre d'arrivée : {delta_v2:.2f} m/s")
print(f"Delta-v total : {delta_v1 + delta_v2:.2f} m/s")

# Calcul de la vitesse héliocentrique à ce point
v_heliocentric = np.sqrt(mu_earth / a_transfer)  # Vitesse à la sortie de la sphère d'influence terrestre
print(f"Vitesse héliocentrique à ce point : {v_heliocentric:.2f} m/s")

# Calcul de l'excentricité
e = (r_jupiter - r_earth) / (r_jupiter + r_earth)
print(f"Excentricité de l'orbite de transfert : {e:.6f}")

# Calcul de l'anomalie excentrique (E)
E = np.arccos((a_transfer - r_earth) / (a_transfer * e))
print(f"Anomalie excentrique à la sortie de la sphère d'influence terrestre : {E:.2f} radians")




# Calcul de la vitesse héliocentrique à la sortie de la sphère d'influence terrestre
v_heliocentric = np.sqrt(2 * delta_v_initial + mu_earth / a1)

# Vitesse radiale par rapport à la Terre (à adapter selon votre mission)
v_T = 10e3  # Exemple : 10 km/s

# Angle entre V_inf_0 et V_T (à adapter selon votre mission)
xi_degrees = 30  # Exemple : 30 degrés
xi_radians = np.deg2rad(xi_degrees)

# Calcul de la vitesse dans le référentiel héliocentrique (V_S0)
V_S0 = np.sqrt(v_heliocentric**2 + v_T**2 + 2 * v_heliocentric * v_T * np.cos(xi_radians))

print(f"Vitesse héliocentrique à la sortie de la sphère d'influence terrestre : {v_heliocentric:.2f} m/s")
print(f"Vitesse radiale par rapport à la Terre : {v_T:.2f} m/s")
print(f"Angle entre V_inf_0 et V_T : {xi_degrees:.2f} degrés")
print(f"Vitesse dans le référentiel héliocentrique (V_S0) : {V_S0:.2f} m/s")
"""