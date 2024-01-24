# Este Script podia necesitar ajuster dependiendo la version de Chimera con la que se trabaje.
#importar modulos necesarios 
import os
import chimera 

#definir la ruta del xtc
xtc_path = "/home/usuario/Escritorio/imagen_ChimeraX/cell/structures/structures"

#definir la ruta del directorio de salida
output_path = "/home/usuario/Desktop/Movies"

#crear el directorio de salida si no existe
if not os.path.exists(output_path):
	os.makedirs(output_path)


#abrir el archivo .xtc en chimera
chimera.openModels.open(xtc_path)

#Definir los parametros de la grabacion
smooth_frames = 5
fps = 30
speed = 1

#Grabarla pelicula
chimera.movie.record(
	smoothFrames=smooth_frames,
	fps = fps,
	 speed=speed,
	 output=output_path
)
