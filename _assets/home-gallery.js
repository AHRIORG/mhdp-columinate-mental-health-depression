(function () {
  const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function initHomeEventsCarousel(root) {
    const track = root.querySelector(".home-events-carousel__track");
    const slides = Array.from(root.querySelectorAll(".home-events-carousel__slide"));
    const items = Array.from(root.querySelectorAll(".home-events-carousel__item"));
    const dots = Array.from(root.querySelectorAll(".home-events-carousel__picture-dot"));
    const prevButton = root.querySelector("[data-carousel-prev]");
    const nextButton = root.querySelector("[data-carousel-next]");
    const slideCount = Math.max(parseInt(root.dataset.slideCount || "1", 10), 1);
    const imageCount = Math.max(parseInt(root.dataset.imageCount || String(items.length), 10), 1);
    const autoplayMs = Math.max(parseInt(root.dataset.autoplayMs || "7000", 10), 3000);

    if (!track || slides.length === 0 || items.length === 0) {
      return;
    }

    const imageToSlide = new Map();
    slides.forEach((slide, slideIndex) => {
      slide.querySelectorAll(".home-events-carousel__item").forEach((item) => {
        const imageIndex = parseInt(item.dataset.imageIndex || "0", 10);
        imageToSlide.set(imageIndex, slideIndex);
      });
    });

    let currentSlide = 0;
    let currentImage = 0;
    let autoplayId = null;

    function setTrackPosition() {
      const offset = currentSlide * (100 / slideCount);
      track.style.transform = `translateX(-${offset}%)`;
    }

    function syncState() {
      items.forEach((item) => {
        const imageIndex = parseInt(item.dataset.imageIndex || "-1", 10);
        item.classList.toggle("is-targeted", imageIndex === currentImage);
      });

      dots.forEach((dot) => {
        const dotImageIndex = parseInt(dot.dataset.targetImage || "-1", 10);
        const isActive = dotImageIndex === currentImage;
        dot.classList.toggle("is-active", isActive);
        if (isActive) {
          dot.setAttribute("aria-current", "true");
        } else {
          dot.removeAttribute("aria-current");
        }
      });

      if (prevButton) {
        prevButton.disabled = slideCount <= 1;
      }

      if (nextButton) {
        nextButton.disabled = slideCount <= 1;
      }
    }

    function goToImage(imageIndex) {
      const normalizedIndex = ((imageIndex % imageCount) + imageCount) % imageCount;
      currentImage = normalizedIndex;
      currentSlide = imageToSlide.get(normalizedIndex) || 0;
      setTrackPosition();
      syncState();
    }

    function goToSlide(slideIndex) {
      const normalizedSlide = ((slideIndex % slideCount) + slideCount) % slideCount;
      currentSlide = normalizedSlide;
      const firstItem = slides[normalizedSlide].querySelector(".home-events-carousel__item");
      if (firstItem) {
        currentImage = parseInt(firstItem.dataset.imageIndex || "0", 10);
      }
      setTrackPosition();
      syncState();
    }

    function startAutoplay() {
      if (reduceMotion || slideCount <= 1) {
        return;
      }

      stopAutoplay();
      autoplayId = window.setInterval(function () {
        goToSlide(currentSlide + 1);
      }, autoplayMs);
    }

    function stopAutoplay() {
      if (autoplayId !== null) {
        window.clearInterval(autoplayId);
        autoplayId = null;
      }
    }

    if (prevButton) {
      prevButton.addEventListener("click", function () {
        goToSlide(currentSlide - 1);
        startAutoplay();
      });
    }

    if (nextButton) {
      nextButton.addEventListener("click", function () {
        goToSlide(currentSlide + 1);
        startAutoplay();
      });
    }

    dots.forEach((dot) => {
      dot.addEventListener("click", function () {
        const targetImage = parseInt(dot.dataset.targetImage || "0", 10);
        goToImage(targetImage);
        startAutoplay();
      });
    });

    root.addEventListener("mouseenter", stopAutoplay);
    root.addEventListener("mouseleave", startAutoplay);
    root.addEventListener("focusin", stopAutoplay);
    root.addEventListener("focusout", function (event) {
      if (!root.contains(event.relatedTarget)) {
        startAutoplay();
      }
    });

    goToImage(0);
    startAutoplay();
  }

  function bootstrap() {
    document.querySelectorAll("[data-home-events-carousel]").forEach(initHomeEventsCarousel);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", bootstrap, { once: true });
  } else {
    bootstrap();
  }
})();
