# CharaPy

CharaPy (Characterization-Python) pretende ser un paquete de Python enfocado en proporcionar las herramientas necesarias para realizar la *Caracterización de la fracción residual (o fracción plus) de fluidos de reservorio* en función a modelos de distribución determinados. 

### Consiste en la distribución y representación de los hidrocarburos que forman parte de la fracción residual como un número conveniente de pseudocomponentes y en la estimación de los parámetros de EoS para cada uno de ellos.

**¿Qué herramientas utiliza?**

1.   Funciones de distribución
2.   Correlación de propiedades
3.   Lumping

Es articulado de manera tal que modelos no previstos en el desarrollo del paquete puedan ser añadidos con facilidad. Así como nuevos desarrollos de correlaciones que permitan caracterizar con diferentes Ecuaciones de Estado (EoS) y criterios varios para la ejecución del lumping. 

# Available properties
------------------------------------------------------
* Distributed component {
                          masa molar
                          densidad
                          porcentaje molar
                          }
* Thermodynamic variables 
(Peng-Robingson y Soave-Redlich-Kwong - In order to include RKPR) {
                                                                   temperatura crítica
                                                                   presión crítica
                                                                   factor acéntrico
                                                                   }
