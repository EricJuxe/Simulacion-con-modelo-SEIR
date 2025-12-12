<div align="center">
Tecnológico Nacional de México
Instituto Tecnológico de Oaxaca
Ingeniería en Sistemas Computacionales

Materia: Simulación
Docente: Martínez Nieto Adelina

Integrantes del equipo

• Juárez Fernández Eric Aarón
• Sebastián Ochoa José Ángel
• Fuentes López Leonardo

Grupo: 5SB
Fecha: 10 de diciembre del año 2025
</div>





## Simulador SEIR de Dengue

Este proyecto es un simulador del dengue usando un modelo matemático llamado SEIR y una interfaz gráfica hecha con Python y Tkinter. Su objetivo es mostrar de manera sencilla cómo puede propagarse el dengue en una población dependiendo de varios factores como contagios, incubación, clima y duración del brote.

## ¿Qué hace este programa?

Permite ingresar datos como:
* Población total
* Personas en incubación
* Personas inmunes
* Días de incubación
* Días en los que una persona contagia
* Duración total de la simulación
* Fuerza del clima en el contagio
* También se pueden cargar parámetros desde archivos Excel o CSV.

## El programa simula día por día cuántas personas hay en cada grupo:

* S: susceptibles (sanas que se pueden enfermar)
* E: expuestas (infectadas pero sin síntomas)
* I: infectadas (las que contagian)
* R: recuperadas o inmunes
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/e9601f66-e498-4f29-a66e-01d931ee7233" />


Muestra dos gráficas principales:
* Evolución de S, E, I y R a lo largo del tiempo
* La tasa de transmisión beta(t), que cambia según el clima
* Muestra datos importantes como:
* Día del pico de contagio
* Mes del pico
* Casos totales estimados
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/7d246364-658c-47e2-b730-95a750900eca" />

## ¿Qué es el modelo SEIR? 

El modelo SEIR divide a la población en cuatro grupos para poder simular cómo avanza una enfermedad:

*S (Susceptibles): personas sanas que sí se pueden enfermar.
*E (Expuestas): personas que ya se contagiaron pero todavía no muestran síntomas.
*I (Infectadas): personas enfermas que sí contagian a otras.
*R (Recuperadas): personas que ya pasaron la enfermedad y ya no se contagian otra vez.

El modelo calcula cómo las personas van pasando de un grupo a otro conforme pasan los días. Así se puede ver si habrá un pico grande, si sube rápido o si baja lento.

## Cómo funciona el clima en el modelo

El contagio del dengue cambia dependiendo del día del año. Para eso se usa una función:

beta(t) = beta0 * (1 + fuerza_estacional * sin(2pit / 365))

En palabras simples:

Cuando el clima favorece a los mosquitos, suben los contagios.
Cuando el clima baja, los contagios también bajan.
Por eso la gráfica de beta(t) se ve como una onda.
Interfaz gráfica

### El programa incluye una ventana con botones y formularios donde puedes:

*Escribir los valores manualmente.
*Cargar un archivo Excel o CSV.
*Ejecutar la simulación.
*Ver las gráficas sin necesidad de programar.
*Cambiar entre la gráfica SEIR y la gráfica beta(t).
*Activar información al pasar el mouse sobre las curvas.
*Cargar parámetros desde Excel o CSV



El programa puede leer parámetros desde un archivo externo. Usa la hoja llamada "Datos" y solo toma la primera fila.
También incluye una opción para generar una plantilla de Excel lista para rellenar.

## Cómo ejecutar el programa

Instalar Python 3.x
Instalar las dependencias:

pip install numpy pandas matplotlib openpyxl mplcursors


Ejecutar el archivo:
python seir_dengue.py


Estructura del proyecto
seir_dengue.py    Archivo principal con el simulador SEIR y la interfaz
