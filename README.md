matlabimagenes
==============
clear all
close all
im_RGB=imread('imagen3.jpg')
im_RGB=im2double(im_RGB);
im_ycbcr=rgb2ycbcr(im_RGB);

%obtengo la crominancia b y r con las que juego para poder detectar los
%ojos de la persona. Para ello represento dichas componentes con valores 
%entre 0 y 255 y posteriormente con la componente cr detecto la cara ya que
%esta la obtiene como el minimo color de la imagen.
im_1=histeq(im_ycbcr(:,:,2));
im_2=histeq(im_ycbcr(:,:,3));
im_21=histeq(im_ycbcr(:,:,1));
im_222=histeq(im_ycbcr(:,:,2));
figure(1), imshow(im_1);
figure(2), imshow(im_2);
%analisis de im_2
cara2=im_2;
for isaias=1:1
    clear longitud
    clear posicion_inicial
    clear posicion_final
    clear cara
    clear im_341
    if isaias~=1
        im_222=contraste2;
        im_21=contraste1;
    end 
    im_34=histeq(cara2);
    im_341=histeq(im_222);
    im_342=histeq(im_21);
    for i=0.2:0.1:0.6
        contraste1=imadjust(im_341,[],[],i);% incremento de contraste con valor i
    end
    for i=0.2:0.1:0.6
        contraste2=imadjust(im_342,[],[],i);% incremento de contraste con valor i
    end
    for i=0.2:0.1:0.6
        contraste=imadjust(im_34,[],[],i);% incremento de contraste con valor i
        figure(3),imshow(contraste),title(['Ajuste de constraste #',int2str(i)],'Color','b');
        disp(i);
    end
    clear im_34
    clear im_2
    clear cara2
    clear acierto
    clear inicio_y
    clear altura
    clear minimos
    
    im_2=contraste;
    clear contraste
    %detecta la cara en la figura        
    for i=floor(length(im_2(:,1))/100):length(im_2(:,1))-floor(length(im_2(:,1))/100)
        longitud(i)=0;
        posicion_inicial(i)=0;
        posicion_final(i)=0;
        for j=floor(length(im_2(1,:))/100):length(im_2(1,:))-floor(length(im_2(1,:))/100);
            if im_2(i,j)>=0.7 && longitud(i)==0
                posicion_inicial(i)=j;
                longitud(i)=longitud(i)+1;
                posicion_final(i)=j;
            elseif im_2(i,j)>=0.7
                posicion_final(i)=j;
            end
        end

        longitud(i)=posicion_final(i)-posicion_inicial(i);

        
        
    end
    longitud3=longitud/2+posicion_inicial;
    longitud2=mean(longitud3);
    varianza=mean(longitud);
    %calculo de altura y posicion inicial del eje y para recortar
    acierto=0;
    for i=floor(length(im_2(:,1))/100):length(im_2(:,1))-floor(length(im_2(:,1))/100)
        inicio_y=0;
        if acierto==0 && posicion_inicial(i)~=0
            acierto=1;
            inicio_y=i;
        end
    end
    %detector de longitud para altura
    for i=floor(length(im_2(:,1))/100):length(im_2(:,1))-floor(length(im_2(:,1))/100)-3
        altura=0;
            if  posicion_inicial(i)~=0
                altura=i;
            end
        
    end
    for i=floor(length(im_2(:,1))/100):length(im_2(:,1))-floor(length(im_2(:,1))/100)
        if posicion_inicial(i)~=0
            minimos=posicion_inicial(i);
        end
    end
    for i=floor(length(im_2(:,1))/100):length(im_2(:,1))-floor(length(im_2(:,1))/100)
        if posicion_inicial(i)~=0 && posicion_inicial(i)<minimos
            minimos=posicion_inicial(i);
        end
    end
    cara2=imcrop(im_2,[(longitud2-varianza/2) inicio_y varianza altura]); 
    cara21=imcrop(contraste1,[(longitud2-varianza/2) inicio_y varianza altura]); 
    cara222=imcrop(contraste2,[(longitud2-varianza/2) inicio_y varianza altura]); 
    
    figure(4), imshow(cara2);
    for i=0.2:0.1:0.6
        contraste1=imadjust(cara21,[],[],i);% incremento de contraste con valor i
    end
        for i=0.2:0.1:0.6
        contraste2=imadjust(cara222,[],[],i);% incremento de contraste con valor i
    end
    for i=0.2:0.1:0.6
        contraste=imadjust(cara2,[],[],i);% incremento de contraste con valor i
        figure(5),imshow(contraste),title(['Ajuste de constraste #',int2str(i)],'Color','b');
        disp(i);
    end
end
im_ycrcb2(:,:,1)=contraste1;
im_ycrcb2(:,:,2)=contraste2;
im_ycrcb2(:,:,3)=contraste;
im_rgbb=ycbcr2rgb(im_ycrcb2);
figure(6),imshow(im_rgbb)
figure(7),imshow(im_ycrcb2(:,:,1))
figure(8),imshow(im_ycrcb2(:,:,2))
figure(9),imshow(im_ycrcb2(:,:,3))
im_grayy=im_ycrcb2(:,:,2);
A=[1 1 1;1 0 1;1 1 1];
filtro=1./9.*A;
filtro_bajo=filter2(filtro,double(im_grayy));
im_edge=edge(filtro_bajo,'canny');
figure(10),imshow(filtro_bajo);
figure(11),imshow(im_edge);
%intentarlo haciendo la media de la longitud +la posicion inicial como
%punto i y sacar el punto j de dicha media y con el centro de la cara
%moverme hacia los lados la distancia media /2 y recortar esa parte solo 1
%vez. despues uso lo del histograma y compruebo si los ojos estan cerrados
%o abiertos
for j=1:3
imag2=bwareaopen(im_edge,j);%N=1,2,5,10,15,20
figure(12),subplot(3,1,j),imshow(imag2),title(['obj removidos menores a: ',int2str(j)],'color','b');
disp(j);
end

for k=1:1:5
     disp(k);
     se=strel('disk',k); %incremento de la figura disk(circulo)con N=1,2,3,4,5,6,7,8,9
     imag3=imclose(imag2,se);
     figure(13),subplot(5,1,k),imshow(imag3),title(['incremento de disk en ',int2str(k)],'color','b');
end
imag4=imfill(imag3,'holes'); 
figure(14),imshow(imag4),title('imagen con obj removidos');
%%
[B,L] = bwboundaries(imag4,'noholes');
    stats = regionprops(L,'all');
    a = regionprops(L, 'area');
    areas = cat(1, a.Area);
    area_max=max(areas);
    indice_area_max=find(areas==area_max);
    boundary_area_max = B{indice_area_max};
    [filas,columnas]=size(boundary_area_max);
    max_x=max(boundary_area_max(:,2));
    min_x=min(boundary_area_max(:,2));
    max_y=max(boundary_area_max(:,1));
    min_y=min(boundary_area_max(:,1));
    imagen_area_max=zeros(max_x,max_y);
    for i=1:filas
        imagen_area_max(boundary_area_max(i,1),boundary_area_max(i,2))=1;
    end
    x= boundary_area_max(:,2);
    y= boundary_area_max(:,1);
    x_media= round((max_x-min_x)/2);
    x_inicial= (min_x + x_media);
    indices = find(x==x_inicial);
    y_correspondientes= y(indices);
    y_inicial = min(y_correspondientes);
    i=1;
        while(x_inicial ~= min_x)
            x_inicial= x_inicial-1;
            x_vector(i)=x_inicial;
            indices = find(x==x_inicial); % bien
            y_correspondientes= y(indices);
            y_inicial = min(y_correspondientes);
            y_vector(i)=y_inicial;
            i=i+1;
        end
    BW = imagen_area_max;
    se = strel('disk',3);
    BW = imclose(BW,se);
    BW = imfill(BW,'holes');
    BW=~BW;
    [filas,columnas] = size(imagen_area_max);
    col = round(columnas/2)-10;
    javi=find(imagen_area_max(:,col));
    row = min(javi);
    figure,
    imshow(im_ycrcb2(:,:,3));
    hold on;
    plot(x_vector,y_vector,'g*','LineWidth',2);
    x = x_vector;
    y = y_vector;
    longitud_x=length(x);
    longitud_y=length(y);
    abc=[x' y' ones(length(x'),1)]\[-(x'.^2+y'.^2)];
    a = abc(1); b = abc(2); c = abc(3);
    xc = -a/2;
    yc = -b/2;
    radio = sqrt((xc^2+yc^2)-c);
    plot(xc,yc,'gx','LineWidth',2);
    theta = 0:0.001:2*pi;
    Xfit = radio*cos(theta) + xc;
    Yfit = radio*sin(theta) + yc;
    plot(Xfit, Yfit);
    centro_y=yc;
    media_y=mean(y_vector);
    if centro_y <= media_y
        disp('IMAGEN 1: OJOS cerrados')
        message = sprintf('Ojos Cerrados');
        text(15,15,message,'Color','y','FontWeight','bold');
    else
        disp('IMAGEN 1: OJOS abiertos')
        message = sprintf('Ojos Abiertos');
        text(15,15,message,'Color','y','FontWeight','bold');
    end

Matlab
