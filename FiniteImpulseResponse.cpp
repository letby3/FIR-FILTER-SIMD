/*
g++ -c FiniteImpulseResponse.cpp -I"C:/SFML-2.5.1-mingw/include"
g++ FiniteImpulseResponse.o -o sfml-app -L"C:/SFML-2.5.1-mingw/lib" -lsfml-graphics -lsfml-window -lsfml-system
*/

#include <SFML/Graphics.hpp>
#include <bits/stdc++.h>

using namespace std;
using namespace sf;

int main()
{
        Image img;
        img.loadFromFile("Test1.jpg");
        Texture texImg;
        texImg.loadFromImage(img);
        Sprite sprImg;
        sprImg.setTexture(texImg);
        sprImg.setPosition(0,0);

        RenderWindow window(sf::VideoMode(675, 450), "SFML works!");
        

        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();
            window.draw(sprImg);
            window.display();
        }

    return 0;
}