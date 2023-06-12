**Estado del sistema** = número de clientes en el sistema. 
**Longitud de la cola** = número de clientes que esperan servicio O estado del sistema menos número de clientes a quienes se les da el servicio.

**N(t)** = número de clientes en el sistema de colas en el tiempo t (t >= 0).

**P_n(t)** = probabilidad de que exactamente n clientes estén en el sistema en el tiempo t, dado el número en el tiempo 0.

**s** = número de servidores (canales de servicio en paralelo) en el sistema de colas.

**Lambda_n** = tasa media de llegadas (número esperado de llegadas por unidad de tiempo) de nuevos clientes cuando hay n clientes en el sistema.

**Miu_n** = tasa media de servicio en todo el sistema (número esperado de clientes que completan su servicio por unidad de tiempo) cuando hay n clientes en el sistema. Nota: Miu_n representa la tasa combinada a la que todos los servidores ocupados (aquellos que están sirviendo a un cliente) logran terminar sus servicios.

**Lambda, Miu, Ro** =  Cuando Lambda_n es constante para toda n, esta constante se denota por Lambda. Cuando la tasa media
de servicio por servidor ocupado es constante para toda n >= 1, esta constante se denota por Miu. (En este caso, Miu_n = s*Miu cuando n >= s, es decir, cuando los s servidores están ocupados.) 

- En estas circunstancias, 1/Lambda y 1/Miu es el tiempo esperado entre llegadas y el tiempo esperado de servicio, respectivamente. Asimismo, Ro = Lambda/(s*Miu) es el factor de utilización de la instalación de servicio, es decir, la fracción esperada de tiempo que los servidores individuales están ocupados.
